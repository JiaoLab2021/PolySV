#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"
import os
import shutil
import logging


# Make dir
class MyMkdir():
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self):
        # If the log handler has not been added yet, add
        if not MyMkdir.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyMkdir.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

    def make(self, dir_path, restart):
        if restart:
            if not os.path.isdir(dir_path):
                os.makedirs(dir_path)
                self.logger.error(f"Make directory: {dir_path}")
            else:
                self.logger.error(f"Directory already exists: {dir_path}, skip ...")
        else:
            if os.path.isdir(dir_path):
                shutil.rmtree(dir_path)
            os.makedirs(dir_path)
            self.logger.error(f"Make directory: {dir_path}")