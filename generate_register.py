#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import re


# # 定义模版文件

# In[2]:


temp_cc =["""#include "register_new.h"
namespace Chroma {
void register_new(bool & foo)
{
  // New measurement registrations go here
  
  ""","""
}
}"""] 
temp_h = ["""#ifndef __register_h__
#define __register_h__
 
// New measurement headers go here

""",
"""
namespace Chroma {
void register_new(bool &);
}

#endif"""]


# # 定位文件

# In[3]:


files = os.listdir("./compile_src/")
files.remove("CMakeLists.txt")
inline_f =[i for i in files if "inline" in i]
heads=[os.path.splitext(i)[0] + ".h" for i in inline_f if ".cc" in os.path.splitext(i)[1]]
inline_f = ["./compile_src/"+i for i in inline_f ]


# # 写入头文件

# In[4]:


w_heas = ["#include"+'"'+i+'"' for i in heads]
w_heas.append("#include"+'"aniso_spectrum_gaugeact2.h"') # 手动添加 purgaug相关的头文件
with open("register_new.h","w") as f:
    f.write(temp_h[0])
    for i in w_heas:
        f.write(i + "\n")
    f.write(temp_h[1])


# # 写入cc文件

# In[5]:


env_list=[]
for i in inline_f:
    with open(i,"r") as f:
        for each_line in f.readlines():
            searchObj = re.search("Inline(.*)Env",each_line)
            if searchObj:
                env_list.append(searchObj.group())
                break
w_env = ["foo &= %s::registerAll();"%i for i in env_list]
w_env.append("foo &= AnisoSpectrumGaugeActEnv2::registerAll();")# 手动添加 purgaug相关的类


# In[6]:


with open("register_new.cc","w") as f:
    f.write(temp_cc[0])
    for i in w_env:
        f.write(i + "\n")
    f.write(temp_cc[1])


# In[ ]:




