# 结构矩阵分析第一次大作业--静力分析大作业
这是第一次大作业的代码，正在测试阶段

```Python
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
 ```
