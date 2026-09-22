"""甲基化对已保存酶切影响的只读分析。

分析规则
--------
* 甲基化位置由请求给出，为 top 链上的 0-based 物理坐标；环形序列上
  坐标对序列长度取模，因此敏感位置跨原点时仍能正确命中。
* 每种酶在请求中给出其识别模式内的"敏感偏移"集合。某个已保存识别
  位点命中的敏感物理坐标为（``k`` = 识别模式长度）：

  - 正向位点（strand "+"）：``start + offset``；
  - 反向位点（strand "-"）：``start + k - 1 - offset``。

  （反向位点在 top 链上的阅读方向为 3'->5'，模式偏移沿 top 链反向
  映射；``matched`` 长度恒等于 ``k``，故无需原模式即可换算。）
* 任一位点的任一敏感物理坐标被甲基化，即判定该位点被**阻断**。
* 位点按物理切口 ``(top_cut, bottom_cut)`` 归组：

  - 切口上仍有至少一个未阻断位点 -> **保留**，并按存活位点重算
    ``enzymes`` / ``sites``（末端类型与突出端不变，物理切口未移动）；
    - 切口上的位点全部被阻断     -> **移除**。
* 片段由保留下来的切口复用 ``digest.build_fragments`` 重建，
  线性/环形拓扑的全部边界规则与完整酶切完全一致。

整个分析只在内存中进行，不读写数据库，原记录保持不变。
"""

from __future__ import annotations

from .digest import build_fragments


def _sensitive_position(site: dict, offset: int, length: int, circular: bool) -> int:
    """把模式内偏移换算为 top 链 0-based 物理坐标（环形坐标取模）。"""
    k = len(site["matched"])
    if site["strand"] == "+":
        pos = site["start"] + offset
    else:
        pos = site["start"] + k - 1 - offset
    return pos % length if circular else pos


def methylation_analysis(
    sequence: str,
    topology: str,
    saved_cuts: list[dict],
    rule_by_enzyme: dict[str, list[int]],
    methylated_positions,
) -> dict:
    """在已保存酶切结果上施加甲基化阻断，返回新的切口分组与重建片段。

    调用方须已完成全部校验（记录存在、坐标不越界、酶名已知、偏移合法）。

    ``rule_by_enzyme`` 为酶名 -> 去重排序后的敏感偏移列表；
    ``methylated_positions`` 为 top 链 0-based 甲基化坐标的可迭代对象。

    返回::

        {"blocked_sites": [...], "retained_cuts": [...],
         "removed_cuts": [...], "fragments": [...]}
    """
    length = len(sequence)
    circular = topology == "circular"
    methylated = set(methylated_positions)

    blocked_sites: list[dict] = []
    retained_cuts: list[dict] = []
    removed_cuts: list[dict] = []

    for cut in saved_cuts:
        surviving_sites: list[dict] = []
        seen_site_keys: set[tuple] = set()  # 同名酶重复提交会产生完全相同的位点，先去重
        for site in cut["sites"]:
            key = (site["enzyme"], site["start"], site["strand"], site["matched"])
            if key in seen_site_keys:
                continue
            seen_site_keys.add(key)

            offsets = rule_by_enzyme.get(site["enzyme"])
            blocking = []
            if offsets:
                for offset in offsets:
                    pos = _sensitive_position(site, offset, length, circular)
                    if pos in methylated:
                        blocking.append(pos)

            if blocking:
                blocked_sites.append(
                    {
                        "enzyme": site["enzyme"],
                        "start": site["start"],
                        "strand": site["strand"],
                        "matched": site["matched"],
                        "top_cut": cut["top_cut"],
                        "bottom_cut": cut["bottom_cut"],
                        "blocked_by_positions": sorted(set(blocking)),
                    }
                )
            else:
                surviving_sites.append(site)

        if surviving_sites:
            # 物理切口仍然存在：末端类型/突出端不变，仅按存活位点重算归属
            retained_cuts.append(
                {
                    "top_cut": cut["top_cut"],
                    "bottom_cut": cut["bottom_cut"],
                    "end_type": cut["end_type"],
                    "overhang": cut["overhang"],
                    "enzymes": sorted({s["enzyme"] for s in surviving_sites}),
                    "sites": sorted(
                        surviving_sites,
                        key=lambda s: (s["enzyme"], s["start"], s["strand"]),
                    ),
                }
            )
        else:
            removed_cuts.append(cut)

    # saved_cuts 本身已按 (top_cut, bottom_cut) 排序，这里显式排序确保结果确定
    retained_cuts.sort(key=lambda c: (c["top_cut"], c["bottom_cut"]))
    removed_cuts.sort(key=lambda c: (c["top_cut"], c["bottom_cut"]))
    blocked_sites.sort(
        key=lambda s: (
            s["top_cut"],
            s["bottom_cut"],
            s["enzyme"],
            s["start"],
            s["strand"],
        )
    )

    fragments = build_fragments(sequence, topology, retained_cuts)

    return {
        "blocked_sites": blocked_sites,
        "retained_cuts": retained_cuts,
        "removed_cuts": removed_cuts,
        "fragments": fragments,
    }
