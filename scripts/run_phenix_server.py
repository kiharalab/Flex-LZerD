#!/usr/bin/env python2

# should be run via phenix.python

import sys
import os
import time

import select

import mmtbx.command_line.geometry_minimization as gm

fifo_path = sys.argv[1]

fifo_flag_file = os.environ["GM_SERVER_FIFO_FLAG_FILE"]

print("fifo_path = %s"%fifo_path)
print("fifo_flag_file = %s"%fifo_flag_file)

while True:
        with open(fifo_path) as fifo:
                while True:
                        # rlist, wlist, xlist, timeout
                        #select.select([fifo], [], [fifo])
                        select.select([fifo], [], [])
                        data = fifo.readline()
                        if data.strip():
                                print("READ LINE: '%s'"%(data.strip()))
                                print("RUNNING...")
                                st = time.time()
                                try:
                                        gm.run(data.split(), log=sys.stdout)
                                except Exception as e:
                                        print(e)
                                os.remove(fifo_flag_file)
                                et = time.time()
                                print("DONE, %.4f sec elapsed"%(et-st))
                                sys.stdout.flush()
                        else:
                                #print("blank line, sleeping")
                                time.sleep(1)
