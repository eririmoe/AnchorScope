from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

# 跟踪是否已经配置过日志
_logging_configured = False
# 跟踪我们自己添加的处理器
_our_handlers: list[logging.Handler] = []


def setup_logging(
    log_file: Optional[Path] = None,
    level: str = "INFO",
    log_format: Optional[str] = None
) -> None:
    """设置日志系统
    
    Args:
        log_file: 日志文件路径，如果为None则只输出到控制台
        level: 日志级别，可选值：DEBUG, INFO, WARNING, ERROR, CRITICAL
        log_format: 日志格式，如果为None则使用默认格式
    """
    global _logging_configured, _our_handlers
    
    # 默认日志格式
    if log_format is None:
        log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    
    # 配置根日志器
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, level.upper()))
    
    # 只移除我们自己添加的处理器，而不是清除所有处理器
    for handler in _our_handlers[:]:
        root_logger.removeHandler(handler)
        _our_handlers.remove(handler)
    
    # 创建格式化器
    formatter = logging.Formatter(log_format, datefmt="%Y-%m-%d %H:%M:%S")
    
    # 控制台处理器
    console_handler = logging.StreamHandler()
    console_handler.setLevel(getattr(logging, level.upper()))
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    _our_handlers.append(console_handler)
    
    # 文件处理器（如果指定了日志文件）
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
        _our_handlers.append(file_handler)
    
    _logging_configured = True


def get_logger(name: str) -> logging.Logger:
    """获取指定名称的日志器
    
    Args:
        name: 日志器名称，通常使用模块名
        
    Returns:
        配置好的日志器实例
    """
    return logging.getLogger(name)
