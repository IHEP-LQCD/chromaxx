#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
import os
import sys
import mylib

files = []
files = ["inline_tests","inline_cluster_dec"]

def GenerateChromaMacros(files):
    """
    根据files自动产生chroma.cc文件中的依赖 
    """
    ChromaMacros = {"Macros_head": "", "Macros_foo": ""}
    tmp = ['#include "%s.h" \n' % i for i in files]
    ChromaMacros["Macros_head"] = "".join(tmp)
    tmp = ["".join([i.capitalize() for i in j.split("_")]) for j in files]
    tmp = ["  foo &= %sEnv::registerAll(); \n"%i for i in tmp]
    ChromaMacros["Macros_foo"] = "".join(tmp)
    return ChromaMacros

def GenerateMakefile(files):
    """
    根据files自动产生Makefile文件
    """
    MakeMacros = {}
    MakeMacros["Macros_head"] = "".join(["\t%s.h \\ \n"%i for i in files ])
    MakeMacros["Macros_obj"] = "".join(["\t%s.o \\ \n"%i for i in files ])
    return MakeMacros



Macros = GenerateChromaMacros(files)
mylib.replace_macros_file("./template_chroma.cc","chroma.cc",Macros)
Macros = GenerateMakefile(files)
mylib.replace_macros_file("template_Makefile","Makefile",Macros)

