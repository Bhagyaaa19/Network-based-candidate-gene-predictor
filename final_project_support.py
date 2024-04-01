#! /usr/bin/env python3
#  -*- coding: utf-8 -*-
#
# Support module generated by PAGE version 8.0
#  in conjunction with Tcl version 8.6
#    Mar 31, 2024 03:05:45 PM IST  platform: Windows NT

import sys
import tkinter as tk
import tkinter.ttk as ttk
from tkinter.constants import *

import final_project

_debug = True  # False to eliminate debug printing from callback functions.


def main(*args):
    '''Main entry point for the application.'''
    global root
    root = tk.Tk()
    root.protocol('WM_DELETE_WINDOW', root.destroy)
    # Creates a toplevel widget.
    global _top1, _w1
    _top1 = root
    _w1 = final_project.Toplevel1(_top1)
    root.mainloop()


# def Show candidate genes(*args):
#     if _debug:
#         print('final_project_support.Show candidate genes')
#         for arg in args:
#             print ('    another arg:', arg)
#         sys.stdout.flush()

if __name__ == '__main__':
    final_project.start_up()
