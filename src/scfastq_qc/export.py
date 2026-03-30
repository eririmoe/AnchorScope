from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


def export_summary_to_csv(summary: dict[str, Any], output_path: Path) -> None:
    """导出摘要到CSV文件
    
    Args:
        summary: 包含QC分析结果的摘要字典
        output_path: CSV文件输出路径
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # 写入基本信息
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Sample Name", summary.get("sample_name", "")])
        writer.writerow(["Total Reads", summary.get("total_reads", 0)])
        writer.writerow(["Total Bases", summary.get("total_bases", 0)])
        writer.writerow(["Mean Read Length", f"{summary.get('mean_read_length', 0):.2f}"])
        writer.writerow(["Median Read Length", summary.get("median_read_length", 0)])
        writer.writerow(["N50", summary.get("n50", 0)])
        writer.writerow(["Min Read Length", summary.get("min_read_length", 0)])
        writer.writerow(["Max Read Length", summary.get("max_read_length", 0)])
        writer.writerow(["Median Read Qscore", f"{summary.get('median_read_qscore', 0):.2f}"])
        writer.writerow(["Long High-Quality Ratio", f"{summary.get('long_high_quality_ratio', 0):.4f}"])
        writer.writerow(["Reversed Read Ratio", f"{summary.get('reversed_read_ratio', 0):.4f}"])
        writer.writerow(["Correct Anchor Order Ratio", f"{summary.get('correct_anchor_order_ratio', 0):.4f}"])
        
        # 写入锚点检测率
        writer.writerow([])
        writer.writerow(["Anchor Detection Ratios"])
        for anchor_name, ratio in summary.get("anchor_detection_ratio", {}).items():
            writer.writerow([anchor_name, f"{ratio:.4f}"])
        
        # 写入结构分类计数
        writer.writerow([])
        writer.writerow(["Structure Classification Counts"])
        for structure, count in summary.get("structure_counts", {}).items():
            writer.writerow([structure, count])
        
        # 写入结构方向性计数
        writer.writerow([])
        writer.writerow(["Structure Orientation Counts"])
        for structure, orientations in summary.get("structure_orientation_counts", {}).items():
            for orientation, count in orientations.items():
                writer.writerow([f"{structure} ({orientation})", count])


def export_anchor_hits_to_csv(read_results: list[dict], output_path: Path) -> None:
    """导出锚点命中详情到CSV文件
    
    Args:
        read_results: 包含每个读取结果的列表
        output_path: CSV文件输出路径
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Read ID", "Read Length", "Read Qscore", "Structure Label", "Is Reversed", 
                         "Anchor Name", "Start", "End", "Mismatches", "Matched Sequence"])
        
        for result in read_results:
            read_id = result.get("read_id", "")
            read_length = result.get("length", 0)
            read_qscore = result.get("read_qscore", 0)
            structure_label = result.get("structure_label", "")
            is_reversed = result.get("is_reversed", False)
            
            for hit in result.get("hits", []):
                writer.writerow([
                    read_id,
                    read_length,
                    f"{read_qscore:.2f}",
                    structure_label,
                    is_reversed,
                    hit.get("anchor_name", ""),
                    hit.get("start", 0),
                    hit.get("end", 0),
                    hit.get("mismatches", 0),
                    hit.get("matched_sequence", "")
                ])


def export_structure_classification_to_csv(summary: dict[str, Any], output_path: Path) -> None:
    """导出结构分类详情到CSV文件
    
    Args:
        summary: 包含QC分析结果的摘要字典
        output_path: CSV文件输出路径
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Structure Class", "Orientation", "Count", "Percentage"])
        
        total_reads = summary.get("total_reads", 0)
        structure_orientation_counts = summary.get("structure_orientation_counts", {})
        
        for structure, orientations in structure_orientation_counts.items():
            for orientation, count in orientations.items():
                percentage = (count / total_reads * 100) if total_reads > 0 else 0
                writer.writerow([structure, orientation, count, f"{percentage:.2f}%"])
