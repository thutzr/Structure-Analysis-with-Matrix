import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import sparse
from scipy.sparse import *

class Joint():
    '''这是节点类'''
    def __init__(self,x,y,GDOF):
        '''
        para::
            x: 节点的x坐标
            y: 节点的y坐标
            GDOF: 节点的整体位移编号，是一个三维向量
        '''
        self.x = x
        self.y = y
        self.gdof = GDOF

class Element():
    '''这是单元类'''
    def __init__(self,JointNo,GlbDOF,Length,CosA,SinA,EI,EA,mass):
        '''
        para::
            JointNo: 两端节点的编号，为二维的向量
            GlbDOF: 单元的整体自由度编号，每个节点有三个编号，故共六个，为六维向量
            Length: 单元的长度
            CosA: 单元与整体x方向夹角的余弦值，在SetElemProp中计算
            SinA: 单元与整体x方向夹角的正弦值，在SetElemProp中计算
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
    def __init__(self,JointNo,LodDOF,LodVal):
        '''
        para::
            JointNo: 节点编号
            LodDOF: 荷载的方向，即与哪个方向的自由度对应
            LodVal: 荷载的大小
        '''
        self.jointNo = JointNo
        self.lodDof = LodDOF
        self.lodVal = LodVal

class ElemLoad():
    '''这是单元均布荷载类'''
    def __init__(self,ElemNo,Indx,Pos,LodVal):
        '''
        para::
            ElemNo: 单元编号
            Indx: 荷载的方向
            Pos: 荷载的位置
            LodVal: 荷载的大小
        '''
        self.elemNo = ElemNo
        self.idx = Indx
        self.pos = Pos
        self.lodVal = LodVal

def SetElemProp(Elem, Joint):
    '''
    para::
        Elem: 包含所有单元的数组
        Joint: 包含所有节点的数组
    return:
        Elem
    '''
    Nelem = len(Elem)
    for i in range(Nelem):
        j1 = Elem[i].jointNo[0]
        j2 = Elem[i].jointNo[1]
        Elem[i].glbdof[:3] = Joint[j1-1].gdof
        Elem[i].glbdof[3:] = Joint[j2-1].gdof
        dx = Joint[j2-1].x - Joint[j1-1].x
        dy = Joint[j2-1].y - Joint[j1-1].y
        Elem[i].len = np.sqrt((dx)**2+(dy)**2)
        Elem[i].cos = dx/Elem[i].len
        Elem[i].sin = dy/Elem[i].len
        Elem[i].mass = 0
    return Elem

def TransMatrix(CosA,SinA):
    '''
    para::
        CosA: 单元与整体坐标系之间的夹角的余弦值
        SinA: 单元与整体坐标系之间的夹角的余弦值
    return:
        ET: 局部坐标系下的单元坐标转换矩阵
    '''
    ET = np.zeros((6,6))
    ET[0,:2] = [CosA, SinA]
    ET[1,:2] = [-SinA, CosA]
    if(ET[0,:].shape[0] > 2):
        ET[2,2] = 1
    if (ET[0,:].shape[0] > 3):
        ET[3:6,3:6] = ET[:3,:3]
    return ET

class Kcol():
    '''这是整体刚度矩阵的一列'''
    def __init__(self,row):
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

    NGlbDOF = len(Kcol[0,:])
    NElem = len(Elem)
    # row1: 每一行的初始行码
    row1 = [NGlbDOF-1]*NGlbDOF
    for i in range(NElem):
        ElocVec = Elem[i].glbdof
        minDOF = 9999999
        for j in range(len(ElocVec)):
            if (ElocVec[j] > 0):
                if (ElocVec[j] < minDOF):
                    minDOF = ElocVec[j]
        for j in range(len(ElocVec)):
            if(ElocVec[j] > 0):
                row1[ElocVec[j]-1] = min([minDOF,row1[ElocVec[j]-1]])
    rowIdx = []
    KVal = []
    indPtr = []
    num = 0
    # 可能有bugs
    for i in range(NGlbDOF):
        rowIdx.append([x for x in range(row1[i],NGlbDOF-row1[NGlbDOF-1-i])])
        temp = NGlbDOF-row1[NGlbDOF-1-i]-row1[i]+1
        indPtr.append([x for x in range(num,num + temp)])
        num = num + temp
        KVal.append(Kcol[row1[i]:temp+row1[i],i])

#    for i in range(NGlbDOF):
#        Kcol[i].row = [0]*Kcol[0,:].length
    return (rowIdx, KVal, indPtr,row1)

def varBandSolv(Disp, Kcol, GLoad, row1):
    NCol = len(Kcol[0,:])
    Diag = np.array([Kcol[i,i] for i in range(NCol)])

   #  for j in range(1,NCol):
        # for k in range(row1[j],j-1):
            # row_1 = max([row1[j],row1[k]])
            # s = np.dot(Kcol[row_1-1:k-1,k],Kcol[row_1-1:k-1,j])
            # Kcol[k,j] = Kcol[k,j]-s
       # # print( Kcol[row1[j]:j-1,j])
        # Kcol[row1[j]-1:j-1,j] = Kcol[row1[j]-1:j-1,j]/Diag[row1[j]-1:j-1]
       # # print( Kcol[row1[j]:j-1,j])
        # s = np.dot(Diag[row1[j]-1:j-1],Kcol[row1[j]-1:j-1,j]**2)
        # Diag[j] = Diag[j] - s
    temp = Kcol.T - np.diag(Kcol.diagonal())
    # print(Kcol.T)
    # print(Kcol.diagonal())
    # print(temp)
    Kcol = Kcol + temp
    # print(Kcol)
    # Kcol += Kcol.T - np.diag(Kcol.diagonal())
    # print(Kcol[:10,:10])
    # print(Kcol)
    # Disp = GLoad
    Disp = np.linalg.solve(Kcol,GLoad)
    # for j in range(1,NCol):
    #     Disp[j] -= np.dot(Kcol[row1[j]-1:j-1,j],Disp[row1[j]-1:j-1])

    # Disp = np.true_divide(Disp,Diag)
    # for j in range(NCol-1,0,-1):
    #     Disp[row1[j]-1:j-1] = Disp[row1[j]-1:j-1] - Disp[j]*Kcol[row1[j]-1:j-1,j]
    return (Kcol, Disp)


def solveDisp(Disp,Elem,Joint,JLoad,ELoad):
    '''
    para::
        Disp: 整体位移向量
        Elem: 包含所有单元的数组
        Joint: 包含所有节点的数组
        JLoad: 节点荷载向量
        ELoad: 单元荷载向量
    return:
        Kcol: 整体刚度矩阵
        Disp: 整体位移向量
    '''
    NGlbDOF = len(Disp)
    GLoad = np.zeros((NGlbDOF))
    Kcol = np.array([[0]*NGlbDOF]*NGlbDOF)
    (rowIdx, KVal, indPtr, row1) = SetMatBand(Kcol, Elem)
    # 得到整体荷载向量
    GLoad = GLoadVec(GLoad, Elem, JLoad, ELoad, Joint)
    # 得到整体刚度矩阵
    Kcol = GStifMat(Kcol, Elem)
    # 得到整体刚度矩阵和位移向量
    (Kcol, Disp) = varBandSolv(Disp, Kcol, GLoad, row1)
    print(Kcol)
    return (Kcol,Disp)

def GStifMat(Kcol, Elem):
    '''
    para::
        Kcol: 整体刚度矩阵
        Elem: 包含所有单元的数组
    return:
        Kcol
    '''
    NElem = len(Elem)
    for i in range(NElem):
        # 计算局部坐标系下的单元刚度矩阵
        EK = EStifMat(Elem[i].len,Elem[i].ei,Elem[i].ea)
        # 计算单元坐标转换矩阵
        ET = TransMatrix(Elem[i].cos,Elem[i].sin)
        # 整体坐标系下的单元刚度矩阵
        EK = np.matmul(np.transpose(ET),np.matmul(EK,ET))
        # 单元定位向量
        ELocVec = Elem[i].glbdof
        for j in range(6):
            JGDOF = ELocVec[j]
            if(JGDOF == 0):
                continue
            # 当节点的位移编码不为0时
            for k in range(len(ELocVec)):
                # 将单元刚度矩阵集成到整体刚度矩阵上
                if((ELocVec[k] > 0) and (ELocVec[k] <= JGDOF)):
                    Kcol[ELocVec[k]-1,JGDOF-1] += EK[k,j]
    # print(Kcol)
    return Kcol

def GLoadVec(GLoad,Elem,JLoad,ELoad,Joint):
    '''
    para::
        GLoad: 整体荷载向量
        Elem: 包含所有单元的数组
        JLoad: 节点荷载向量
        ELoad: 单元荷载向量
        Joint: 包含所有节点的向量
    return:
        GLoad: 整体荷载向量
    '''
    NJLoad = len(JLoad)
    for i in range(NJLoad):
        JGDOF = Joint[JLoad[i].jointNo].gdof[JLoad[i].lodDof]
        GLoad[JGDOF] = GLoad[JGDOF]+JLoad[i].lodVal
    NELoad = len(ELoad)
    for i in range(NELoad):
        ET = TransMatrix(Elem[ELoad[i].elemNo].cos,Elem[ELoad[i].elemNo].sin)
        F0 = EFixendF(ELoad[i].idx,ELoad[i].pos,ELoad[i].lodVal,Elem[ELoad[i].elemNo])
        F0 = np.matmul(np.transpose(ET),F0)
        ELocVec = np.array(Elem[ELoad[i].elemNo].glbdof)
        nonZeroIdx = np.where(ELocVec>0)[0]
        GLoad[ELocVec[nonZeroIdx]-1] += F0[nonZeroIdx]

    return GLoad

def EStifMat(ELen,EI,EA):
    '''
    para::
        Elem: 包含所有单元的数组
        EI: 单元的抗弯刚度
        EA: 单元的抗拉刚度
    return:
        EK: 单元刚度矩阵
    '''
    EK = np.array([[0]*6]*6)
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

def EFixendF(Indx,a,q,Elem):
    '''
    para::
        Indx: 力的方向
        a: 力的位置
        q: 力的大小
    return:
        F0: 单元力
    '''
    l = Elem.len
    EI = Elem.ei
    EA = Elem.ea
    F0 = np.zeros((6))
    if(Indx == 1):
        F0[1] = -q*(a*l)/2*(2-2*(a*l)**2/(l**2)+(a*l)**3/(l**3))
        F0[2] = -q*(a*l)**2 /12 *(6-8*(a*l)/l+3*(a*l)**2/(l**2))
        F0[4] = -q*(a*l)**3/(l**2)/2*(2-(a*l)/l)
        F0[5] = q*(a*l)**3/l/12*(4-3*(a*l)/l)
    elif(Indx == 2):
        F0[1] = -q*(l-(a*l))**2/(l**2)*(1+2*(a*l)/l)
        F0[2] = -q*(a*l)*(l-(a*l))**2/(l**2)
        F0[4] = -q*(a*l)**2/(l**2)*(1-2*(l-(a*l))/l)
        F0[5] = q*(a*l)**2*(l-(a*l))/(l**2)
    elif(Indx == 3):
        F0[0] = EA*q/l
        F0[3] = -EA*q/l
    elif(Indx == 4):
        F0[1] = 12*EI*q/(l**3)
        F0[2] = 6*EI*q/(l**2)
        F0[4] = -12*EI*q/(l**3)
        F0[5] = 6*EI*q/(l**2)

    for i in range(F0.shape[0]):
        F0[i] = -F0[i]

    return F0

def ElemDisp(Disp,Elem):
    '''
    para::
        Disp: 整体位移向量
        Elem: 一个单元类
    return:
        EDisp: 单元位移向量
    '''
    ET = TransMatrix(Elem.cos,Elem.sin)
    EDisp = np.zeros((6))
    for i in range(len(Elem.glbdof)):
        if(Elem.glbdof[i] > 0):
            EDisp[i] += Disp[Elem.glbdof[i]-1]
    EDisp = np.matmul(np.transpose(ET),EDisp)

    return EDisp

def ElemForce(ie, Disp,Elem,ELoad):
    '''
    para::
        ie: 单元的编号
        Disp: 整体位移向量
        Elem: 一个单元
        ELoad: 单元荷载
    return:
        EForce: 单元荷载向量
    '''
    ET = TransMatrix(Elem.cos,Elem.sin)
    EK = EStifMat(Elem.len,Elem.ei,Elem.ea)
    EDisp = ElemDisp(Disp,Elem)
    # 总的单元荷载向量
    F = [0]*6
    for i in range(len(ELoad)):
        if(ELoad[i].elemNo == ie):
            # 由一个单元荷载引起的单元荷载向量
            F0 = EFixendF(ELoad[i].idx,ELoad[i].pos,ELoad[i].lodVal,Elem)
            F = F + F0

    EForce = np.matmul(EK,np.matmul(ET,EDisp))-F
    EForce = EForce*np.array([-1,1,1,1,1,1])

    return EForce


def InputData():
    with open('SM90.IPT','r') as f:
        PropType = f.readline().strip('\n')
        if(int(PropType) != 1):
            return
        [NElem,NJoint,NGlbDOF,NJLoad,NELoad] = f.readline().strip('\n').split(',')
        NElem = int(float(NElem))
        NJoint = int(float(NJoint))
        NGlbDOF = int(NGlbDOF)
        NJLoad = int(NJLoad)
        NELoad = int(NELoad)
        jointList = []
        Elem = []
        JLoad = []
        ELoad = []

        for i in range(NJoint):
            temp = f.readline().strip('\n').split(',')
            x = float(temp[0])
            y = float(temp[1])
            gdof = [int(s) for s in temp[2:]]
            joint = Joint(x,y,gdof)
            jointList.append(joint)

        for i in range(NElem):
            [jointN1,jointN2,EA,EI] = f.readline().strip('\n').split(',')
            EA = int(float(EA.replace('d','e')))
            # 全部初始化为0
            elem = Element([int(float(jointN1)),int(float(jointN2))],[],0,0,0,int(float(EI)),EA,0)
            Elem.append(elem)
        if(NJLoad > 0):
            for j in range(NJLoad):
                [jointN,lodDof,lodVal] = f.readline().split(',')
                jointLod = JointLoad(int(jointN),int(lodDof),int(float(lodVal)))
                JLoad.append(jointLod)
        if(NELoad > 0):
            for j in range(NELoad):
                [elemNo, Indx, Pos, LodVal] = f.readline().strip('\n').split(',')
                elemLoad = ElemLoad(int(elemNo),int(Indx),float(Pos),int(float(LodVal)))
                ELoad.append(elemLoad)

    return (NElem, NJoint,NGlbDOF,NJLoad,NELoad,jointList,Elem,JLoad,ELoad)


def OutputResult(NElem,Disp,Elem,ELoad):
    with open('SMCAI90.out','w') as f:
        f.write('10 0\n')
        EDisp = [0]*6
        for i in range(NElem):
           EDisp = ElemDisp(Disp, Elem[i])
           for j in range(len(ELoad)):
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
           for k in range(6):
               if(k < 5):
                   f.write(str(EDisp[k]))
                   f.write(' ')
               else:
                   f.write(str(EDisp[k]))
                   f.write('\n')
        for i in range(NElem):
            EForce = ElemForce(i, Disp, Elem[i], ELoad)
            for k in range(len(EForce)):
                if(k < len(EForce)-1):
                    f.write(str(EForce[k]))
                    f.write(' ')
                else:
                    f.write(str(EForce[k]))
                    f.write('\n')
    return

def main():
    (NElem, NJoint,NGlbDOF,NJLoad,NELoad,
            jointList,Elem,JLoad,ELoad) = InputData()
    Elem = SetElemProp(Elem, jointList)

    Disp = [0]*NGlbDOF
    (Kcol,Disp) = solveDisp(Disp,Elem,jointList,JLoad,ELoad)
    OutputResult(NElem,Disp,Elem,ELoad)
    return

main()
