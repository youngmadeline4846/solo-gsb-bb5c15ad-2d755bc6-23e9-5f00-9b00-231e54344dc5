"""甲基化影响分析接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


# 非回文 4bp 模式、平端切于中间：GGCC -> cut(2,2)；
# 反向位点命中模式 revcomp(GGCC)=GGCC 是回文，故测试反链使用 GTCA。
# GTCA 非回文（revcomp=TGAC），反向命中时 start+4-cut_top 等换算便于手算。
BLUNT4 = {"name": "Blunt4", "pattern": "GGCC", "cut_top": 2, "cut_bottom": 2}
REV4 = {"name": "Rev4", "pattern": "GTCA", "cut_top": 2, "cut_bottom": 2}
ECO = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}


def _create(client, sequence, topology="linear", enzymes=None):
    payload = {
        "sequence": sequence,
        "topology": topology,
        "enzymes": enzymes if enzymes is not None else [ECO],
    }
    resp = client.post("/digests", json=payload)
    assert resp.status_code == 201
    return resp.json()


def _analyze(client, digest_id, positions, rules):
    return client.post(
        f"/digests/{digest_id}/methylation-analysis",
        json={"methylated_positions": positions, "enzyme_rules": rules},
    )


def _rule(name, *offsets):
    return {"enzyme": name, "sensitive_offsets": list(offsets)}


def test_forward_site_blocked_cut_removed_and_fragments_rebuilt(client):
    # AAAA GAATTC TTTT：EcoRI 位点 start=4，offset 1 -> 物理位置 5
    body = _create(client, "AAAAGAATTCTTTT")
    resp = _analyze(client, body["id"], [5], [_rule("EcoRI", 1)])
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["id"] == body["id"]
    assert data["sequence_length"] == 14
    assert data["methylated_positions"] == [5]

    assert len(data["blocked_sites"]) == 1
    blocked = data["blocked_sites"][0]
    assert blocked == {
        "enzyme": "EcoRI",
        "start": 4,
        "strand": "+",
        "matched": "GAATTC",
        "top_cut": 5,
        "bottom_cut": 9,
        "blocked_by_positions": [5],
    }
    assert data["retained_cuts"] == []
    assert len(data["removed_cuts"]) == 1
    assert (data["removed_cuts"][0]["top_cut"],
            data["removed_cuts"][0]["bottom_cut"]) == (5, 9)
    assert [f["length"] for f in data["fragments"]] == [14]


def test_non_sensitive_offset_does_not_block(client):
    # offset 0 映射到位置 4；甲基化位置 5 不在敏感集合中 -> 不阻断
    body = _create(client, "AAAAGAATTCTTTT")
    resp = _analyze(client, body["id"], [5], [_rule("EcoRI", 0)])
    assert resp.status_code == 200
    data = resp.json()
    assert data["blocked_sites"] == []
    assert len(data["retained_cuts"]) == 1
    assert data["removed_cuts"] == []
    assert [f["length"] for f in data["fragments"]] == [5, 9]


def test_methylation_outside_any_site_is_noop(client):
    body = _create(client, "AAAAGAATTCTTTT")
    resp = _analyze(client, body["id"], [0, 13], [_rule("EcoRI", 1)])
    assert resp.status_code == 200
    data = resp.json()
    assert data["blocked_sites"] == []
    assert len(data["retained_cuts"]) == 1
    assert data["methylated_positions"] == [0, 13]


def test_shared_physical_cut_retained_with_surviving_site(client):
    # 两个酶各一个正向位点（GGCC start=5 与 GTCA start=13），均为平端，
    # 产生两个不同物理切口 (7,7) 与 (15,15)；仅甲基化阻断 GGCC 位点
    # （offset 1 -> 位置 6），验证保留/移除切口分组与片段重建。
    body = _create(client, "AAAARGGCCTTTTGTCAA", enzymes=[BLUNT4, REV4])
    assert [c["top_cut"] for c in body["cuts"]] == [7, 15]

    resp = _analyze(client, body["id"], [6], [_rule("Blunt4", 1)])
    assert resp.status_code == 200
    data = resp.json()
    assert [s["enzyme"] for s in data["blocked_sites"]] == ["Blunt4"]
    retained = [(c["top_cut"], c["enzymes"]) for c in data["retained_cuts"]]
    assert retained == [(15, ["Rev4"])]
    assert [c["top_cut"] for c in data["removed_cuts"]] == [7]
    assert [f["length"] for f in data["fragments"]] == [15, 3]


def test_shared_physical_cut_merged_from_two_sites_retained(client):
    # 同一序列上两条不同模式的位点都产生物理切口 (8,8)：
    # GGCC start=4 与一条反向 GTCA 位点，换算后切口同为 (8,8)。
    # 序列: AAAA GGCC TT? 需要反向 TGAC 起点 start 满足 start+4-2=8 -> start=6
    # 即 top 链 6..9 = TGAC（revcomp(GTCA)）。令序列 4..9 = GGTGAC 与 GGCC 冲突，
    # 故改为两个不同 start 的等价构造：直接放 GGCC 在 4，TGAC 在 6 需要
    # seq[6:10]=TGAC；seq[4:8]=GGCC -> seq[4..7]=GGCC，故 seq[6..9]=CCTG，
    # 无法同时满足。退而求其次：回文酶也以同一切口几何切于 start=4 之后 4bp。
    #
    # 可行构造：两种酶识别同一 4bp 区段 GGCC（start=5），各自一条正向位点，
    # 均切 (7,7) -> 一个物理切口、两个 site。
    eco_alias = {"name": "Blunt4Alias", "pattern": "GGCC", "cut_top": 2, "cut_bottom": 2}
    body = _create(
        client, "AAAARGGCCTTTT", enzymes=[BLUNT4, eco_alias]
    )
    assert len(body["cuts"]) == 1
    cut = body["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (7, 7)
    assert cut["enzymes"] == ["Blunt4", "Blunt4Alias"]
    assert len(cut["sites"]) == 2

    # 仅阻断 Blunt4 的位点：物理切口由 Blunt4Alias 位点保住，归属被重算
    resp = _analyze(client, body["id"], [6], [_rule("Blunt4", 1)])
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["enzyme"] == "Blunt4"
    assert len(data["retained_cuts"]) == 1
    retained = data["retained_cuts"][0]
    assert retained["top_cut"] == 7
    assert retained["enzymes"] == ["Blunt4Alias"]
    assert [s["enzyme"] for s in retained["sites"]] == ["Blunt4Alias"]
    assert data["removed_cuts"] == []
    assert [f["length"] for f in data["fragments"]] == [7, 6]

    # 两条位点都被阻断 -> 物理切口移除
    resp2 = _analyze(
        client,
        body["id"],
        [6],
        [_rule("Blunt4", 1), _rule("Blunt4Alias", 1)],
    )
    data2 = resp2.json()
    assert len(data2["blocked_sites"]) == 2
    assert data2["retained_cuts"] == []
    assert len(data2["removed_cuts"]) == 1
    assert [f["length"] for f in data2["fragments"]] == [13]


def test_reverse_strand_offset_mapping(client):
    # AAATGACTTTT：top 链 3..6=TGAC = revcomp(GTCA)，Rev4 反向位点
    # start=3，平端切口 start+4-2=5 -> (5,5)。
    # 模式偏移到 top 链物理位置：start + 4 - 1 - offset = 6 - offset。
    # offset 0（模式 5' 端碱基，top 链位点 3' 端）-> 物理位置 6
    body = _create(client, "AAATGACTTTT", enzymes=[REV4])
    assert body["cuts"][0]["sites"][0]["strand"] == "-"
    assert (body["cuts"][0]["top_cut"], body["cuts"][0]["bottom_cut"]) == (5, 5)

    resp = _analyze(client, body["id"], [6], [_rule("Rev4", 0)])
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    blocked = data["blocked_sites"][0]
    assert blocked["strand"] == "-"
    assert blocked["blocked_by_positions"] == [6]
    assert len(data["removed_cuts"]) == 1
    assert [f["length"] for f in data["fragments"]] == [11]

    # 甲基化位置 3（offset 3 映射到 6-3=3）对仅敏感于 offset 0 的规则无效
    resp2 = _analyze(client, body["id"], [3], [_rule("Rev4", 0)])
    assert resp2.json()["blocked_sites"] == []

    # offset 3 -> 6 - 3 = 3，此时阻断
    resp3 = _analyze(client, body["id"], [3], [_rule("Rev4", 3)])
    assert len(resp3.json()["blocked_sites"]) == 1


def test_multiple_sensitive_offsets_any_one_blocks(client):
    # EcoRI 位点 start=4：offset 0/5 分别映射到 4/9
    body = _create(client, "AAAAGAATTCTTTT")
    resp = _analyze(client, body["id"], [9], [_rule("EcoRI", 0, 5)])
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["blocked_by_positions"] == [9]

    # 两个敏感位置同时甲基化时一并回显（排序）
    resp2 = _analyze(client, body["id"], [9, 4], [_rule("EcoRI", 5, 0)])
    assert resp2.json()["blocked_sites"][0]["blocked_by_positions"] == [4, 9]


def test_circular_forward_position_modulo(client):
    # 环形 ATTCGGGGGA（长度 10）：GAATTC 跨原点位点 start=8，
    # offset 1 -> 9；offset 2 -> 10 ≡ 0（跨原点）。
    body = _create(client, "ATTCGGGGGA", topology="circular")
    assert body["cuts"][0]["sites"][0]["start"] == 8
    assert body["cuts"][0]["top_cut"] == 9  # 环形单切口开环

    resp = _analyze(client, body["id"], [0], [_rule("EcoRI", 2)])
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["blocked_by_positions"] == [0]
    assert len(data["removed_cuts"]) == 1
    # 环形无切口：保留整条环状分子
    assert [(f["start"], f["length"]) for f in data["fragments"]] == [(0, 10)]

    # 位置 9 对 offset 2 不敏感（它是 offset 1）-> 不阻断
    resp2 = _analyze(client, body["id"], [9], [_rule("EcoRI", 2)])
    assert resp2.json()["blocked_sites"] == []
    assert len(resp2.json()["retained_cuts"]) == 1


def test_circular_reverse_position_modulo(client):
    # 环形 CGGGGGGGTGA（长度 10）：反向 GTCA 位点 TGAC 跨原点，
    # top[7:10]=TGA、top[0]=C -> start=7，平端切口 7+4-2=9 -> (9,9)。
    # 敏感物理坐标：offset 0 -> start+3 = 10 ≡ 0（跨原点取模）；
    # offset 3 -> start+0 = 7。
    body = _create(client, "CGGGGGGTGA", topology="circular", enzymes=[REV4])
    assert body["sequence_length"] == 10
    sites = body["cuts"][0]["sites"]
    assert len(sites) == 1 and sites[0]["strand"] == "-"
    assert sites[0]["start"] == 7
    assert body["cuts"][0]["top_cut"] == 9

    # offset 0 映射到 10 % 10 = 0（敏感物理坐标跨原点取模）
    resp = _analyze(client, body["id"], [0], [_rule("Rev4", 0)])
    assert resp.status_code == 200
    assert len(resp.json()["blocked_sites"]) == 1
    assert resp.json()["blocked_sites"][0]["blocked_by_positions"] == [0]

    # offset 3 -> 7，不跨原点
    resp2 = _analyze(client, body["id"], [7], [_rule("Rev4", 3)])
    assert len(resp2.json()["blocked_sites"]) == 1

    # 位置 1 对 offset 0（映射到 0）不敏感
    resp3 = _analyze(client, body["id"], [1], [_rule("Rev4", 0)])
    assert resp3.json()["blocked_sites"] == []
    assert len(resp3.json()["retained_cuts"]) == 1

    # 阻断唯一环形切口 -> 保留整条环状分子
    resp4 = _analyze(client, body["id"], [0], [_rule("Rev4", 0)])
    assert [f["length"] for f in resp4.json()["fragments"]] == [10]


def test_multiple_sites_partial_blocking_linear(client):
    # 两个 EcoRI 位点：GAATTC at 0 and 10
    seq = "GAATTCAAAAGAATTC"
    body = _create(client, seq)
    assert [c["top_cut"] for c in body["cuts"]] == [1, 11]
    # 阻断第一个位点（offset 1 -> 位置 1）
    resp = _analyze(client, body["id"], [1], [_rule("EcoRI", 1)])
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["start"] == 0
    assert [c["top_cut"] for c in data["retained_cuts"]] == [11]
    assert [c["top_cut"] for c in data["removed_cuts"]] == [1]
    assert [f["length"] for f in data["fragments"]] == [11, 5]


def test_blocked_sites_and_cuts_are_stably_sorted(client):
    seq = "GAATTCAAAAGAATTC"
    body = _create(client, seq)
    # 同时阻断两个位点，响应中甲基化坐标与 blocked_sites 均稳定排序；
    # 与请求传入顺序无关
    r1 = _analyze(client, body["id"], [11, 1], [_rule("EcoRI", 1)]).json()
    r2 = _analyze(client, body["id"], [1, 11], [_rule("EcoRI", 1)]).json()
    assert r1 == r2
    assert r1["methylated_positions"] == [1, 11]
    assert [s["start"] for s in r1["blocked_sites"]] == [0, 10]
    keys = [(s["top_cut"], s["bottom_cut"]) for s in r1["blocked_sites"]]
    assert keys == sorted(keys)
    assert r1["retained_cuts"] == []


def test_rule_for_enzyme_without_sites_is_allowed(client):
    # EcoRV 在序列中没有位点；给出它的规则不报错，也不阻断任何东西
    ecov = {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}
    body = _create(client, "AAAAGAATTCTTTT", enzymes=[ECO, ecov])
    resp = _analyze(
        client,
        body["id"],
        [5],
        [_rule("EcoRI", 1), _rule("EcoRV", 0)],
    )
    assert resp.status_code == 200
    assert [s["enzyme"] for s in resp.json()["blocked_sites"]] == ["EcoRI"]


def test_missing_record_returns_404(client):
    resp = _analyze(client, 999, [0], [_rule("EcoRI", 0)])
    assert resp.status_code == 404


def test_methylated_position_out_of_bounds_returns_422(client):
    body = _create(client, "AAAAGAATTCTTTT")  # 长度 14
    for bad in (14, 100):
        resp = _analyze(client, body["id"], [bad], [_rule("EcoRI", 1)])
        assert resp.status_code == 422, resp.text
        assert "out of bounds" in resp.text


def test_negative_position_returns_422(client):
    body = _create(client, "AAAAGAATTCTTTT")
    resp = client.post(
        f"/digests/{body['id']}/methylation-analysis",
        json={"methylated_positions": [-1], "enzyme_rules": [_rule("EcoRI", 0)]},
    )
    assert resp.status_code == 422


def test_unknown_enzyme_name_returns_422(client):
    body = _create(client, "AAAAGAATTCTTTT")
    resp = _analyze(client, body["id"], [5], [_rule("NopEI", 0)])
    assert resp.status_code == 422
    assert "unknown enzyme" in resp.text


def test_offset_out_of_range_returns_422(client):
    # EcoRI 模式长度 6：合法偏移 0..5
    body = _create(client, "AAAAGAATTCTTTT")
    for bad in (6, 7):
        resp = _analyze(client, body["id"], [5], [_rule("EcoRI", bad)])
        assert resp.status_code == 422, resp.text
        assert "out of range" in resp.text


def test_duplicate_positions_rules_offsets_return_422(client):
    body = _create(client, "AAAAGAATTCTTTT")

    r1 = client.post(
        f"/digests/{body['id']}/methylation-analysis",
        json={"methylated_positions": [5, 5], "enzyme_rules": [_rule("EcoRI", 1)]},
    )
    assert r1.status_code == 422
    assert "duplicate methylated" in r1.text

    r2 = client.post(
        f"/digests/{body['id']}/methylation-analysis",
        json={
            "methylated_positions": [5],
            "enzyme_rules": [_rule("EcoRI", 1), _rule("EcoRI", 2)],
        },
    )
    assert r2.status_code == 422
    assert "duplicate enzyme rules" in r2.text

    r3 = client.post(
        f"/digests/{body['id']}/methylation-analysis",
        json={"methylated_positions": [5], "enzyme_rules": [_rule("EcoRI", 1, 1)]},
    )
    assert r3.status_code == 422
    assert "duplicate sensitive offsets" in r3.text


def test_empty_lists_and_missing_fields_return_422(client):
    body = _create(client, "AAAAGAATTCTTTT")
    base = f"/digests/{body['id']}/methylation-analysis"

    assert client.post(
        base, json={"methylated_positions": [], "enzyme_rules": [_rule("EcoRI", 1)]}
    ).status_code == 422
    assert client.post(
        base, json={"methylated_positions": [5], "enzyme_rules": []}
    ).status_code == 422
    assert client.post(base, json={}).status_code == 422


def test_multiple_validation_problems_aggregated(client):
    body = _create(client, "AAAAGAATTCTTTT")
    resp = client.post(
        f"/digests/{body['id']}/methylation-analysis",
        json={
            "methylated_positions": [99],
            "enzyme_rules": [_rule("NopEI", 0), _rule("EcoRI", 6)],
        },
    )
    assert resp.status_code == 422
    assert "out of bounds" in resp.text
    assert "unknown enzyme" in resp.text
    assert "out of range" in resp.text


def test_blank_enzyme_name_returns_422(client):
    body = _create(client, "AAAAGAATTCTTTT")
    resp = client.post(
        f"/digests/{body['id']}/methylation-analysis",
        json={"methylated_positions": [5],
              "enzyme_rules": [{"enzyme": "  ", "sensitive_offsets": [0]}]},
    )
    assert resp.status_code == 422


def test_methylation_analysis_does_not_persist(client):
    body = _create(client, "AAAAGAATTCTTTT")
    _analyze(client, body["id"], [5], [_rule("EcoRI", 1)])
    # 原记录内容不变，且没有新增任何提交
    stored = client.get(f"/digests/{body['id']}").json()
    assert stored == body
    assert len(client.get("/digests").json()) == 1
