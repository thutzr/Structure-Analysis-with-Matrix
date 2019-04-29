import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import sparse
from scipy.sparse import *

class Joint():
    '''这是节点类'''
    def __init__(x,y,GDOF):
        '''
        para::
            x: 节点的x坐标
            y: 节点的y坐标
            GDOF: 节点的整体位移编号
        '''
        self.x = x
        self.y = y
        self.gdof = GDOF

class Element():
    '''这是单元类'''
    def __init__(JointNo,GlbDOF,Length,CosA,SinA,EI,EA,mass):
        '''
        para::
            JointNo: 两端节点的编号，为二维的向量
            GlbDOF: 单元的整体自由度编号，每个节点有三个编号，故共六个，为六维向量
            Length: 单元的长度
            CosA: 单元与整体x方向夹角的余弦值
            SinA: 单元与整体x方向夹角的正弦值
            EI: 单元抗弯刚度
            EA: 单元的抗拉刚度
            mass: 单元的质量
        '''
        self.jointNo = JointNo
        self.glbdof = GlbDOF
        self.len = Length
        self.cos = CosA
        self.sin = SinA
        self.ei = EI
        self.ea = EA
        self.m = mass

class JointLoad():
    '''这是节点荷载类'''
    def __init__(JointNo,LodDOF,LodVal):
        '''
        para::
            JointNo: 节点编号
            LodDOF: 荷载的方向，即与哪个方向的自由度对应
            LodVal: 荷载的大小
        '''
        self.jointNo = JointNo
        self.lodDof = lodDOF
        self.lodVal = LodVal

class ElemLoad():
    '''这是单元荷载类'''
    def __init__(ElemNo,Indx,Pos,LodVal):
        '''
        para::
            ElemNo: 单元编号
            Indx: 荷载的方向？
            Pos: 荷载的位置
            LodVal: 荷载的大小
        '''
        self.elemNo = ElemNo
        self.idx = Indx
        self.pos = Pos
        self.lodVal = LodVal

def SetElemProp(Elem, Joint):
    Nelem = Elem.length
    for i in range(Nelem):
        j1 = Elem[i].jointNo[0]
        j2 = Elem[i].jointNo[1]
        Elem[i].glbdof[0:3] = Joint[j1].gdof
        Elem[i].glbdof[3:6] = Joint[j2].gdof
        dx = Joint[j2].x - Joint[j1].x
        dy = Joint[j2].y - Joint[j1].y
        Elem[i].len = sqrt(dx**2+dy**2)
        Elem[i].cos = dx/Elem[i].len
        Elem[i].sin = dy/Elem[i].len
        Elem[i].mass = 0
    return Elem
def TransMatrix(ET,CosA,SinA):
    '''
    para::
        ET:
        CosA: 单元与整体坐标系之间的夹角的余弦值
        SinA: 单元与整体坐标系之间的夹角的余弦值
    '''
    ET[0,:2] = [CosA, SinA]
    ET[1,:2] = [-SinA, CosA]
    if(ET[0,:].length > 2):
        ET[2,2] = 1
    if (ET[0,:].length > 3):
        ET[3:6,3:6] = ET[:3,:3]
    return ET

class Kcol():
    '''这是整体刚度矩阵的一列'''
    def __init__(row):
        '''
        para::
            row:
        '''
        self.row = row

# 采用CSC稀疏矩阵存储方式
def SetMatBand(Kcol,Elem):
    '''
    para::
        Kcol: 整体刚度矩阵，与上述定义的Kcol不同
        Elem: 向量
    '''

    NGlbDOF = Kcol[0,:].length
    NElem = Elem.length
    # row1: 每一行的初始行码
    row1 = [NGlbDOF-1]*NGlbDOF

    for i in range(NElem):
        ElocVec = Elem[i].GlbDOF
        minDOF = 9999999
        for j in range(ElocVec.length):
            if (ElocVec[j] > 0):
                if (ElocVec[j] < minDOF):
                    minDOF = ElocVec[j]
        for j in range(ElocVec.length):
            if(ElocVec[j] > 0):
                row1[ElocVdc[j]-1] = min([minDOF,row1[ElocVec[j]-1]])
    rowIdx = []
    KVal = []
    indPtr = []
    num = 0
    # 可能有bugs
    for i in range(NGlbDOF):
        rowIdx.append([x for x in range(row1[i],NGlbDOF-row1[NGlDOF-1-i]))
        temp = NGlbDOF-row1[NGlbDOF-1-i]-row1[i]+1
        indPtr.append([x for x in range(num,num + temp)])
        num = num + temp
        KVal.append(Kcol[row1[i]:temp+row1[i],i])

#    for i in range(NGlbDOF):
#        Kcol[i].row = [0]*Kcol[0,:].length
    return (rowIdx, KVal, indPtr,row1)

def varBandSolv(Disp, Kcol, GLoad,row1):
    NCol = Kcol[0,:].length
    Diag = np.array([Kcol[i].row[i] for i in range(NCol)])
    for j in range(1,NCol):
        for k in range(row1[j],j-1):
            row_1 = max([row1[j],row1[k]])
            s = np.dot(Kcol[k].row[row_1:k-1],Kcol[j].row[row_1:k-1])
            Kcol[j].row[k] = Kcol[j].row[k]-s
        Kcol[j].row[row1[j]:j-1] = Kcol[j].row[row1[j]:j-1]/Diag[row[j]:j-1]
        s = np.dot(Diag[row1[j]:j-1],Kcol[j].row[row1[j]:j-1]**2)
        Diag[j] = Diag[j] - s
    Disp = GLoad
    for j in range(1,NCol):
        Disp[j] = Disp[j] - np.dot(Kcol[j].row[row1[j]:j-1],Disp[row1[j]:j-1])
    Disp = Disp/Diag
    for j in range(NCol-1,0,-1):
        Disp[row1[j]:j-1] = Disp[row1[j]:j-1] - Disp[j]*Kcol[j]/row[row1[j]:j-1]

    return (Kcol, Disp)


def solveDisp(Disp,Elem,Joint,JLoad,ELoad):
    Disp = [0]*Disp.length
    NGlbDOF = Disp.length
    (rowIdx, KVal, indPtr,row1) = SetMatBand(Kcol,Elem)
    GLoad = GLoadVec(GLoad,Elem,JLoad,Eload,Joint)
    Kcol = GStifMat(Kcol,Elem)
    (Kcol, Disp) = varBandSolv(Disp,Kcol,GLoad)
    return Disp

def GStifMat(Kcol, Elem):
    NElem = Elem.length
    for i in range(NElem):
        EStifMat(EK,Elem[i].len,Elem[i].ei,Elem[i].ea)
        ET = TransMatrix(ET,Elem[i].cos,Elem[i].sin)
        EK = np.matmul(np.transpose(ET),np.matmul(EK,ET))
        ELocVec = Elem[i].glbdof
        for j in range(6):
            JGDOF = ELocVec(j)
            if(JGDOF == 0):
                continue
            for k in range(ELocVec.length):
                if((ELocVec[k] > 0) and (ELocVec[k] <= JGDOF)):
                    Kcol[JGDOF].row[ELocVec[k]] = Kcol[JGDOF].row[ELocVec[k]] + EK[:,j]
    return Kcol

def GLoadVec(GLoad,Elem,JLoad,ELoad,Joint):
    GLoad = [0]*GLoad.length
    NJLoad = Jload.length
    for i in range(NJLoad):
        JGDOF = Joint[JLoad[i].jointNo].GDOF[JLoad[i].lodDof]
        GLoad[JGDOF] = GLoad[JGDOF]+JLoad[i].lodVal

    NEload = ELoad.length
    for i in range(NELoad):
        ET = TransMatrix(ET,Elem[ELoad[i].elemNo].cos,Elem[ELoad[i].elemNo].sin ))
        EFixendF(F0,ELoad[i].idx,ELoad[i].pos,ELoad[i].lodVal,Elem[ELoad[i].ElemNo])
        F0 = np.matmul(np.transpose(ET),F0)
        ELocVec = Elem(ELoad[i].elemNo).glbdof
        nonZeroIdx = np.where(ELocVec>0)
        GLoad[ELocVec[nonZeroIdx]] += F0

    return GLoad

def EStifMat(EK,Elen,EI,EA):
    EAL = EA/ELen
    EIL1 = EI/ELen
    EIL2 = EI/(ELen**2)
    EIL3 = EI/(ELen**3)
    EK[0,:] = [EAL,0,0,-EAL,0,0]
    EK[1,:] = [0,12*EIL3,6*EIL2,0,-12*EIL3,6*EIL2]
    EK[2,:] = [0,6*EIL2,4*EIL1,0,-6*EIL2,2*EIL1]
    EK[3,:] = [-EAL,0,0,EAL,0,0]
    EK[4,:] = [0,-12*EIL3,-6*EIL2,0,12*EIL3,-6*EIL2]
    EK[5,:] = [0,6*EIL2,2*EIL1,0,-6*EIL2,4*EIL1]
    return EK

def EFixendF(F0,Indx,a,q,Elem):
    F0 = [0]*F0.length
    l = Elem.len
    EI = Elem.ei
    EA = Elem.ea
    if(Indx == 1):
        F0[1] = -q*(a*l)/2*(2-2*(a*l)**2/l)
        F0[2] = -q*(a*l)**2 /12 *(6-8*(a*l)/l+3*(a*l)**2/(l**2))
        F0[4] = -q*(a*l)**3/(l**2)/2*(2-(a*l)/l)
        F0[5] = q*(a*l)**3/l/12*(4-3*(a*l)/l)
    else if(Indx == 2):
        F0[1] = -q*(l-(a*l))**2/(l**2)*(1+2*(a*l)/l)
        F0[2] = -q*(a*l)*(l-(a*l))**2/(l**2)
        F0[4] = -q*(a*l)**2/(l**2)*(1-2*(l-(a*l))/l)
        F0[5] = q*(a*l)**2*(l-(a*l))/(l**2)
    else if(Indx == 3):
        F0[0] = EA*q/l
        F0[3] = -EA*q/l
    else if(Indx == 4):
        F0[1] = 12*EI*q/(l**3)
        F0[2] = 6*EI*q/(l**2)
        F0[4] = -12*EI*q/(l**3)
        F0[5] = 6*EI*q/(l**2)

    F0 = -F0
    return F0

def ElemDisp(EDisp,Disp,Elem):
    EDisp = [0]*EDisp.length
    ET = TransMatrix(ET,Elem.cos,Elem.sin)
    EDisp = np.matmul(np.transpose(ET),EDisp)
    nonZeroIdx = np.where(Elem.glbdof)
    EDisp[nonZeroIdx] = EDisp[nonZeroIdx] + Disp[Elem.glbdof[nonZeroIdx]]

    return EDisp

def ElemForce(ie, Disp,Elem,ELoad):
    ET = TransMatrix(ET,Elem.cos,Elem.sin)
    EK = EStifMat(EK,Elem.len,Elem.ei,Elem.ea)
    EDisp = ElemDisp(EDisp,Disp,Elem)
    F = [0]*6
    for i in range(ELoad.length):
        if(ELoad[i].elemNo == ie):
            F0 = EFixendF(F0,ELoad[i].idx,ELoad[i].pos,ELoad[i].lodVal,Elem)
            F = F + F0

    EForce = np.matmul(EK,np.matmul(E,EDisp))-F
    EForce = EForce*np.array([-1,1,1,1,1,1])
    return EForce

def InputData():
    with open('SM90.IPT','r') as f:
        PropType = f.readline()
        if(PropType != 1):
            return
        [NElem,NJoint,NGlbDOF,NJLoad,NELoad] = f.readline().split(',')
        NElem = int(NElem)
        NJoint = int(NJoint)
        NGlbDOF = int(NGlbDOF)
        NJLoad = int(NJLoad)
        NELoad = int(NELoad)
        jointList = []
        Elem = []
        JLoad = []
        ELoad = []

        for i in range(NJoint):
            [x,y,gdof] = f.readline().split(',')
            joint = Joint(float(x),float(y),int(gdof))
            jointList.append(joint)

        for i in range(NElem):
            [jointN1,jointN2,EA,EI] = f.readline().split(',')
            EA = int(float(EA.replace('d','e')))
            # 全部初始化为0
            elem = Element([int(jointN1),int(jointN2)],[],0,0,0,int(EI),EA,0)
            Elem.append(elem)
        if(NJLoad > 0):
            [jointN,lodDof,lodVal] = f.readline().split(',')
            jointLod = JointLoad(int(jointN),int(lodDof),int(float(lodVal)))
            JLoad.append(jointLod)
        if(NELoad > 0):
            [elemNo, Indx, Pos, LodVal] = f.readline().split(',')
            elemLoad = ElemLoad(int(elemNo),int(Indx),float(Pos),int(float(LodVal)))
            ELoad.append(elemLoad)

    return (NElem, NJoint,NGlbDOF,NJLoad,NELoad,jointList,Elem,JLoad,ELoad)


def OutputResult(NElem,EDisp,Disp,Elem,ELoad,EForce):
    with open('SMCAI90.out','w') as f:
        f.write('10 0\n')
        for i in range(NElem):
           EDisp = ElemDisp(ElemDisp, Disp, Elem[i])
           for j in range(ELoad.length):
               if((ELoad[j].idx == 3) and (ELoad[j].elemNo == i+1)):
                   if(ELoad[j].pos == 0):
                       EDisp[0] = EDisp[0] + ELoad[j].lodVal
                    else:
                        EDisp[3] = EDisp[3] + ELoad[j].lodVal
                if((ELoad[j].idx==4) and (ELoad[j].elemNo==i+1)):
                    if(ELoad[j].pos == 0):
                        EDisp[1] = EDisp[1] + ELoad[j].lodVal
                    else:
                        EDisp[4] = EDisp[4] + ELoad[j].lodVal
            f.write(EDisp)
        for i in range(NElem):
            EForce = ElemForce(i, Disp, Elem[i], ELoad)
            f.write(EForce)
    return

def main():
    (NElem, NJoint,NGlbDOF,NJLoad,NELoad,jointList,Elem,JLoad,ELoad) = InputData()
    Elem = SetElemProp(Elem, jointList)
    Disp = [0]*NGlbDOF
    Disp = solveDisp(Disp,Elem,jointList,JLoad,ELoad)
    OutputResult(NElem,EDisp,Disp,Elem,ELoad,EForce)
    return

main()
