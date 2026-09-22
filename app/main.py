"""FastAPI 应用入口。

启动时（lifespan）自动完成数据库初始化；接口：
    POST /digests                          提交序列与酶列表，执行酶切并保存，返回完整结果
    GET  /digests                          列出全部提交的摘要
    GET  /digests/{id}                     按 id 查询保存的完整酶切结果
    POST /digests/{id}/variant-analysis    基于已保存记录施加单碱基替换，重跑酶切并对比差异（不落库）
    POST /digests/{id}/enzyme-panel-plan   枚举候选酶组合，寻找达到目标片段数的最少酶集合（只读，不落库）
    POST /digests/{id}/partial-digest-distribution
                                           为每个已保存切口指定独立切割概率，枚举成功/失败
                                           状态并返回片段长度多重集的概率分布（只读，不落库）
    GET  /health                           健康检查
"""

from __future__ import annotations

import os
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException

from .db import Database
from .digest import run_digest
from .distribution import MAX_PARTIAL_CUTS, partial_digest_distribution
from .panel import plan_enzyme_panel
from .schemas import (
    DigestRequest,
    DigestResponse,
    DigestSummary,
    EnzymeIn,
    EnzymePanelPlanRequest,
    EnzymePanelPlanResponse,
    PartialDigestDistributionResponse,
    PartialDigestRequest,
    VariantAnalysisRequest,
    VariantAnalysisResponse,
)
from .variant import apply_substitutions, diff_cuts, diff_fragments, variant_summary

db: Database


@asynccontextmanager
async def lifespan(app: FastAPI):
    global db
    db = Database(os.environ.get("DATABASE_PATH", "data/digest.db"))
    yield


app = FastAPI(
    title="DNA Restriction Digest Simulator",
    version="1.0.0",
    description="DNA 限制性酶切模拟：支持线性/环形序列、IUPAC 简并碱基、"
                "双链切口偏移、黏端/平端分析与片段生成。",
    lifespan=lifespan,
)


def _build_response(record: dict) -> dict:
    return {
        "id": record["id"],
        "topology": record["topology"],
        "sequence": record["sequence"],
        "sequence_length": len(record["sequence"]),
        "enzymes": record["enzymes"],
        "cuts": record["result"]["cuts"],
        "fragments": record["result"]["fragments"],
        "created_at": record["created_at"],
    }


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/digests", status_code=201, response_model=DigestResponse)
def create_digest(request: DigestRequest) -> dict:
    # Pydantic 已在此前完成全部校验（含切口偏移合法性），
    # 非法请求在写库之前即以 422 拒绝，保证整次保存被拒绝。
    result = run_digest(request.sequence, request.topology, request.enzymes)
    submission_id, created_at = db.save_submission(
        request.sequence, request.topology, request.enzymes, result
    )
    return _build_response(
        {
            "id": submission_id,
            "sequence": request.sequence,
            "topology": request.topology,
            "created_at": created_at,
            "enzymes": [e.model_dump() for e in request.enzymes],
            "result": result,
        }
    )


@app.get("/digests", response_model=list[DigestSummary])
def list_digests() -> list[dict]:
    return db.list_submissions()


@app.get("/digests/{submission_id}", response_model=DigestResponse)
def get_digest(submission_id: int) -> dict:
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")
    return _build_response(record)


@app.post(
    "/digests/{submission_id}/variant-analysis",
    response_model=VariantAnalysisResponse,
)
def analyze_variants(submission_id: int, request: VariantAnalysisRequest) -> dict:
    # 先确认记录存在：不存在返回 404，且不泄露任何校验细节
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")

    sequence = record["sequence"]
    length = len(sequence)
    # 位置越界依赖记录的序列长度，只能在取到记录后校验 -> 422
    out_of_bounds = sorted(
        v.position for v in request.variants if v.position >= length
    )
    if out_of_bounds:
        raise HTTPException(
            status_code=422,
            detail=(
                f"variant positions out of bounds for sequence length {length}: "
                + ", ".join(str(p) for p in out_of_bounds)
            ),
        )

    # 仅在内存中施加替换并重跑酶切，原记录与数据库均不改动
    enzymes = [EnzymeIn(**e) for e in record["enzymes"]]
    variant_sequence = apply_substitutions(sequence, request.variants)
    before = record["result"]
    after = run_digest(variant_sequence, record["topology"], enzymes)

    return {
        "id": record["id"],
        "sequence_length": length,
        "variant_sequence": variant_sequence,
        "variants": variant_summary(sequence, request.variants),
        "cuts": diff_cuts(before["cuts"], after["cuts"]),
        "fragments": diff_fragments(before["fragments"], after["fragments"]),
    }


@app.post(
    "/digests/{submission_id}/enzyme-panel-plan",
    response_model=EnzymePanelPlanResponse,
)
def plan_panel(submission_id: int, request: EnzymePanelPlanRequest) -> dict:
    # 先确认记录存在：不存在返回 404，且不泄露任何校验细节
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")

    # 仅在内存中枚举候选组合并复用酶切计算，原记录与数据库均不改动
    plan = plan_enzyme_panel(
        record["sequence"],
        record["topology"],
        request.candidates,
        request.target_fragment_count,
    )
    if plan is None:
        raise HTTPException(
            status_code=422,
            detail=(
                f"no combination of the {len(request.candidates)} candidate "
                f"enzyme(s) yields exactly {request.target_fragment_count} "
                "fragment(s)"
            ),
        )

    return {
        "id": record["id"],
        "target_fragment_count": request.target_fragment_count,
        "selected_enzymes": [e.model_dump() for e in plan["enzymes"]],
        "fragment_count": plan["fragment_count"],
        "cuts": plan["cuts"],
    }


@app.post(
    "/digests/{submission_id}/partial-digest-distribution",
    response_model=PartialDigestDistributionResponse,
)
def partial_digest(submission_id: int, request: PartialDigestRequest) -> dict:
    # 先确认记录存在：不存在返回 404，且不泄露任何校验细节
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")

    saved_cuts = record["result"]["cuts"]

    # 枚举规模有界：超过 12 个物理切口无法在保证耗时的前提下全枚举 -> 422
    if len(saved_cuts) > MAX_PARTIAL_CUTS:
        raise HTTPException(
            status_code=422,
            detail=(
                f"submission {submission_id} has {len(saved_cuts)} saved cuts; "
                f"partial-digest enumeration supports at most {MAX_PARTIAL_CUTS}"
            ),
        )

    saved_keys = {(c["top_cut"], c["bottom_cut"]) for c in saved_cuts}
    requested = {(item.top_cut, item.bottom_cut): item.probability for item in request.cuts}

    # 每个已保存物理切口必须恰好出现一次：缺失或多余均 -> 422
    missing = sorted(saved_keys - set(requested))
    extra = sorted(set(requested) - saved_keys)
    if missing or extra:
        problems = []
        if missing:
            problems.append(
                "missing probabilities for cuts: "
                + ", ".join(f"({t}, {b})" for t, b in missing)
            )
        if extra:
            problems.append(
                "probabilities given for unknown cuts: "
                + ", ".join(f"({t}, {b})" for t, b in extra)
            )
        raise HTTPException(status_code=422, detail="; ".join(problems))

    # 概率按已保存切口的顺序（(top_cut, bottom_cut) 升序）排列，
    # 与请求条目的传入顺序无关，保证枚举与结果确定
    probabilities = [
        requested[(c["top_cut"], c["bottom_cut"])] for c in saved_cuts
    ]

    # 仅在内存中枚举成功/失败状态并重建片段，原记录与数据库均不改动
    distribution = partial_digest_distribution(
        record["sequence"], record["topology"], saved_cuts, probabilities
    )

    return {
        "id": record["id"],
        "topology": record["topology"],
        "sequence_length": len(record["sequence"]),
        "cut_count": len(saved_cuts),
        "outcome_count": len(distribution["outcomes"]),
        "total_probability": distribution["total_probability"],
        "outcomes": distribution["outcomes"],
    }
