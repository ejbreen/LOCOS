# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 18:48:54 2017

@author: Evan
"""

import threading
import time
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

lock = threading.Lock()

a, b, c, d, e, f = range(2), range(4), range(3), range(4), range(5), range(2)

thingies = [a, b, c, d, e, f]

def printEverything(thingy):
    time.sleep(len(thingy))
    with lock:
         print thingy
    
pool = ThreadPool(4)

pool.map(printEverything, thingies)

pool.close()

pool.join()



threads = 6
length = 144
seg = int(round(length/threads,0))


segs = []
for i in reversed(list(range(threads))):
    print i
    temp = list(range(int(round(length-((i+1)*seg),0)), 
                      int(round(length-(i*seg),0))))
    segs.append(temp)
    print temp
    print segs
