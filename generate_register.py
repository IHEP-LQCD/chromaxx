#!/usr/bin/env python
# coding: utf-8

# In[14]:


import os
import re


# # 定义模版文件

# In[61]:


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

# In[16]:


files = os.listdir("./compile_src/")
files.remove("CMakeLists.txt")
inline_f =[i for i in files if "inline" in i]
heads=[os.path.splitext(i)[0] + ".h" for i in inline_f if ".cc" in os.path.splitext(i)[1]]
inline_f = ["./compile_src/"+i for i in inline_f ]


# # 写入头文件

# In[42]:


w_heas = ["#include"+'"'+i+'"' for i in heads]
with open("register_new.h","w") as f:
    f.write(temp_h[0])
    for i in w_heas:
        f.write(i + "\n")
    f.write(temp_h[1])


# # 写入cc文件

# In[48]:


env_list=[]
for i in inline_f:
    with open(i,"r") as f:
        for each_line in f.readlines():
            searchObj = re.search("Inline(.*)Env",each_line)
            if searchObj:
                env_list.append(searchObj.group())
                break
w_env = ["foo &= %s::registerAll();"%i for i in env_list]


# In[63]:


with open("register_new.cc","w") as f:
    f.write(temp_cc[0])
    for i in w_env:
        f.write(i + "\n")
    f.write(temp_cc[1])


# In[ ]:




