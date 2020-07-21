# -*- coding: utf-8 -*-
"""
Created on Tue Dec  14 15:17 2019

@author: wpiper
Time practice
"""
import serial
import time
import datetime

def IRon(ser):
    ser.write(b'A')
    print('IR on!')
def IRoff(ser):
    ser.write(b'B')
    print('IR off!')
def FIBERon(ser):
    ser.write(b'C')
    print('FIBER on!')
def FIBERoff(ser):
    ser.write(b'D')
    print('FIBER off!')
def AUDIOon(ser):
    ser.write(b'E')
    print('AUDIO on!')
def AUDIOoff(ser):
    ser.write(b'F')
    print('AUDIO off!')
def SHOCKERon(ser):
    ser.write(b'G')
    print('SHOCKER on!')
def SHOCKERoff(ser):
    ser.write(b'H')
    print('SHOCKER off!')
    
def StartProt(ser):
    IRon(ser)
    FIBERon(ser)
    print('PROTOCOL START')
    print(datetime.datetime.now() -tzero)
def StartProtEnd(ser):
    IRoff(ser)
    FIBERoff(ser)
def CS(ser):
    IRon(ser)
    FIBERon(ser)
    AUDIOon(ser)
    print('CS on')
    print(datetime.datetime.now() -tzero)
def CSint(ser):
    IRoff(ser)
    FIBERoff(ser)
def US(ser):
    IRon(ser)
    FIBERon(ser)
    SHOCKERon(ser)
    print('US on')
    print(datetime.datetime.now() -tzero)
def USend(ser):
    IRoff(ser)
    FIBERoff(ser)
    SHOCKERoff(ser)
    AUDIOoff(ser)
    print('CS and US off')
    print(datetime.datetime.now() -tzero)
def CSend(ser):
    IRoff(ser)
    FIBERoff(ser)
    AUDIOoff(ser)
    print('CS off')
    print(datetime.datetime.now() -tzero)
def UUSend(ser):
    IRoff(ser)
    FIBERoff(ser)
    SHOCKERoff(ser)
    print('US off')
    print(datetime.datetime.now() -tzero)
def EndProtStart(ser):
    IRon(ser)
    FIBERon(ser)
def EndProtEnd(ser):
    IRoff(ser)
    FIBERoff(ser)
    print(datetime.datetime.now() -tzero)
    



#tfunclist = [StartProt,StartProtEnd,CS,CSint,US,USend,EndProtStart,EndProtEnd]

funclist = [StartProt,StartProtEnd]

startendbuffertime = 4
acctime = 300 - startendbuffertime
cooldowntime = 150 - startendbuffertime

tlist=[1]
tlist.append(tlist[-1]+startendbuffertime)

#first ITItimes is acclimation time
#ITI times must be integers!!!
ITItimes = iter([acctime,200,150,190,180,120,150,135])#]200,150,190,180,150])#,210,160,140,170,220]) #middle 150 is before U-US
#300 CS1, 530 CS2, 710 CS3, 930 CS4, 1140 CS5, 1320 UUS, 1470
for i in range(8):
    if i < 5:
        funclist.append(CS)
        tlist.append(tlist[-1]+next(ITItimes))
        funclist.append(CSend)
        tlist.append(tlist[-1]+30) #CS length
    if i >=5:
        funclist.append(US)
        tlist.append(tlist[-1]+next(ITItimes))
        funclist.append(UUSend)
        tlist.append(tlist[-1]+5) #UUS length
funclist.append(EndProtStart)
funclist.append(EndProtEnd)
tlist.append(tlist[-1]+cooldowntime)
tlist.append(tlist[-1]+startendbuffertime)

print(len(funclist))
print(funclist)
print(len(tlist))
flist=[]
for t in tlist:
    flist.append(t-1)
print(flist)
tfunclist = funclist
endtime = tlist[-1]+5
print(endtime)

for i,j in zip(tlist, funclist): print(i-1), print(j), print('\n')

with serial.Serial('COM4', 9600, timeout=1) as ser:
    
    t0 = time.time()  
    timechecked = False
    timer = 0
    tzero = datetime.datetime.now()
    time.sleep(2)
    ser.write(b'B')
    ser.write(b'D')
    ser.write(b'F')
    ser.write(b'H')
    t1 = time.time()
    delta_t0 = t1 - t0
    while delta_t0 < endtime:
        t1 = time.time()
        delta_t0 = t1 - t0
        if timechecked == False:
            t2 = time.time()
            timechecked = True
        delta_t1 = t1 - t2
        if delta_t1 > 0.999:
            timer += 1
            #print(delta_t1)
            #print(time.time())
            #print(datetime.datetime.now()-tzero)
            print(timer)
            timechecked = False
            if timer in tlist:
                a = tlist.index(timer)
                tfunclist[a](ser)
    ser.write(b'B')
    ser.write(b'D')
    ser.write(b'F')
    ser.write(b'H')
        
print('COMPLETE')

"""

with serial.Serial('COM5', 9600, timeout=1) as ser:   
    time.sleep(2)
    ser.write(b'G')
    ser.write(b'A')
    ser.write(b'C')
    time.sleep(5)
    ser.write(b'H')
    ser.write(b'B')
    ser.write(b'D')
"""