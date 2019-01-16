#!/usr/bin/env python
# encoding: utf-8
'''
@author: FENG
@contact: WEI_Lingfeng@163.com
@file: DC2M.py
@time: 2019/1/11 15:33
@desc:
'''

import re
import numpy as np
import cv2
from matplotlib import pyplot as plt
import pandas as pd

class RMC_ANA:
    #count=0
    def __init__(self,name,scale):
        self.name=name#计算名称
        self.scale=scale+2#栅格行列数，相等情况下

    '''
    def count(self):
        print('Total RMC_ANA %d' % (RMC_ANA.count))
    '''

    def read_tally(self):
        flname=self.name+'.Tally'
        fp = open(flname, 'r')  # 读取文件
        text = fp.read()
        return text
        fp.close()

    def POW(self):
        L=self.scale
        text=self.read_tally()
        Np=L*L#总珊元数
        '燃耗计算的功率分布数据接口'
        #编辑正则表达式
        pattp=['']*Np
        for i in range(Np):
            pattp[i] = ">."+str(i+1)+".>.0.*(\d\.\d*E.\d*).*\d\..*\n"
        pattbu="Total.Burnup..MWD/KgHM.*\d\.\d{3}"
        bu=re.findall(pattbu,text)

        #匹配正则表达式
        lenbu=len(bu)
        busteps=int(len(bu)/2)
        POW_cell=['']*lenbu#每个燃耗步celli的功率的集合
        POW=[POW_cell]*Np
        for i in range(Np):
            POW[i]=re.findall(pattp[i],text)


        #格式化处理
        fp = open(self.name + '.POW', 'w')
        for i in range(lenbu):
            fp.write(bu[i]+'\n')
            fp.write('\n')
            for j in range(L):
                str1=''
                    #print(POW_L[j][i])
                    #if j>0:  print(POW_L[j-1][i])

                for k in range(L):
                    str1+=POW[j*L+k][i]+'\t'
                    #print(str1)
                fp.write(str1+'\n')
            fp.write('\n')
        fp.close()

    def NPOW(self):
        'normalized power'
        L = self.scale
        text = self.read_tally()
        Np = L * L  # 总珊元数
        # 编辑正则表达式
        pattp = [''] * Np
        for i in range(Np):
            pattp[i] = ">." + str(i + 1) + ".>.0.*(\d\.\d*E.\d*).*\d\..*\n"
        pattbu = "Total.Burnup..MWD/KgHM.*(\d\.\d{3})"
        bu = re.findall(pattbu, text)

        # 匹配正则表达式
        if bu is None:
            # 字符串编写
            patt1 = ">.1.>.*>.*(\d\.\d*E.\d*).*\d\..*\n"
            patt2 = ">.2.>.*>.*(\d\.\d*E.\d*).*\d\..*\n"
            patt4 = ">.4.>.*>.*(\d\.\d*E.\d*).*\d\..*\n"
            patt5 = ">.5.>.*>.*(\d\.\d*E.\d*).*\d\..*\n"
            # "(\d\.\dE-\d*).*\d\.\d*\n"
            m1 = re.findall(patt1, text)
            m2 = re.findall(patt2, text)
            m4 = re.findall(patt4, text)
            m5 = re.findall(patt5, text)
            # 转成浮点类型
            M1 = [0] * len(m1)
            M2 = [0] * len(m1)
            M4 = [0] * len(m1)
            M5 = [0] * len(m1)
            for i in range(len(m1)):
                M1[i] = float(m1[i])
                M2[i] = float(m2[i])
                M4[i] = float(m4[i])
                M5[i] = float(m5[i])
            Ptotal = sum(M1) + sum(M2) + sum(M4) + sum(M5)
            Pave = Ptotal / (len(m1) - 25) / 4

            # 换算成归一化功率
            for i in range(len(m1)):
                m1[i] = '{:.5f}'.format(M1[i] / Pave)
                m2[i] = '{:.5f}'.format(M2[i] / Pave)
                m4[i] = '{:.5f}'.format(M4[i] / Pave)
                m5[i] = '{:.5f}'.format(M5[i] / Pave)

            print((sum(M1)) / Pave)

            # 拼接12组的同行
            m12 = [(' ')] * 17
            # print(type(m12[1]))
            for i in range(17):
                for j in range(17):
                    m12[i] += m1[i * 17 + j] + ' '
                for k in range(17):
                    m12[i] += m2[i * 17 + k] + ' '

            # 拼接45组的同行
            m45 = [(' ')] * 17
            for i in range(17):
                for j in range(17):
                    m45[i] += m4[i * 17 + j] + ' '
                for k in range(17):
                    m45[i] += m5[i * 17 + k] + ' '

            # 输出文件
            fp = open(flname + '.txt', 'w')
            for i in range(len(m12)):
                fp.write(m12[i] + '\n')
            for i in range(len(m45)):
                fp.write(m45[i] + '\n')
            fp.close()#

        else:
            lenbu = len(bu)
            POW_cell = [''] * lenbu  # 每个燃耗步cell_i的功率的集合
            POW = [POW_cell] * Np

            POWff_cell = [''] * lenbu  # 格式化的
            POWff = [POW_cell] * Np
            for i in range(Np):
                POW[i] = re.findall(pattp[i], text)

            POWf=np.ones((Np,lenbu),dtype=float)
            for i in range(lenbu):
                for j in range(Np):
                    POWf[j][i]=float(POW[j][i])
            Psum=POWf.sum(0)

            POWS=np.zeros([(L+1)*lenbu,L],dtype=float)
            for i in range(lenbu):
                POWS[(L+1)*i][0]=float(bu[i])
                for j in range(L):
                    for k in range(L):
                        POWS[(L+1)*i+j+1][k]=POWf[L*j+k][i]/Psum[i]*312
            index = []
            for i in range(lenbu):
                index += ['burnup']
                for j in range(L):
                    index += [str(i+1)+'0'+str(j + 1)]
            df=pd.DataFrame(POWS,index)
            df.to_csv(self.name+'.csv',encoding='utf-8')
            return df
            '''
            #Normalization step
            for i in range(lenbu):
                for j in range(Np):
                    POWf[j][i] = POWf[j][i]/Psum[i]*312
            POWff_cell = [0] * lenbu
            POWff = [POWff_cell] * Np
            POWfl = POWf.tolist()
            fp = open(self.name + '.NPOW', 'w')
            for i in range(lenbu):
                fp.write(bu[i] + '\n')
                fp.write('\n')
                for j in range(L):

                    for k in range(L):
                        POWff[L*j+k][i]= '{:.3f}'.format(POWfl[L*j+k][i])
                        fp.write(str(POWff[L*j+k][i]) + '\t')
                    fp.write('\n')
                fp.write('\n')
            fp.close()
            '''

    def plt_POW(self,burnup_step):

        import os

        bustep=burnup_step*2
        ximg=1000 #图像长和宽
        yimg=1000

        D=int(yimg*0.8)  #组件对边长
        L=D*2/np.sqrt(3) #组件对角线长
        d=2/np.sqrt(3)/(self.scale-2+1/3)*D
        l=2/np.sqrt(3)*d
        xmid=int(ximg/2)
        ymid=int(yimg/2)

        "组件边界"
        Marray=[[xmid+L/2,ymid],[xmid+L/4,ymid+D/2],[xmid-L/4,ymid+D/2],[xmid-L/2,ymid],[xmid-L/4,ymid-D/2],[xmid+L/4,ymid-D/2]]
        img = np.zeros((ximg, yimg, 3), np.uint8)
        Mpts = np.array(Marray,np.int32)
        #mpts = mpts.reshape((-1, 1, 2))
        cv2.polylines(img, [Mpts], True, (0, 255, 255),2)

        "珊元绘图"
        #第一个珊元
        L2=int(self.scale/2)*d
        fxmid=xmid-np.sqrt(3)*L2*np.sqrt(3)/2#中点
        fymid=ymid-np.sqrt(3)*L2/2
        marray = [[fxmid + d / 2, fymid - l / 4], [fxmid + d / 2, fymid + l / 4], [fxmid, fymid + l / 2],
                  [fxmid - d / 2, fymid + l / 4], [fxmid - d / 2, fymid - l / 4], [fxmid, fymid - l / 2]]#第一个珊元各角点坐标
        vector1=np.array([[d,0]]*6)#横向平移向量
        vector2=np.array([[d/2,l*3/4]]*6)#斜向平移向量
        if os.path.exists(self.name + '.csv') is True:
            print('归一化功率分布数据csv文件已读取成功。')
            df = pd.read_csv(self.name + '.csv',header=0,index_col=0)
        else:print('文件未存在，请先执行RMC_ANA.NPOW()方法。')


        lochead=str(bustep)+'01'#切片索引开头
        locend=str(bustep)+'023'#切片索引结尾
        plt_data=df.loc[lochead:locend]#切片

        plt_data_np = np.array(plt_data)#转成np数据类型，以便计算
        plt_data=pd.DataFrame(plt_data_np)
        #由于np小数位数不能完美指定（例如：要求0.000，但在np中只能达到0.0的效果），
        #而pandas数据格式能够，因此下面使用pandas格式定义格式后再转成text字符串数组，
        #这里，plt_data重新赋值的原因为：重新定义索引成0-22，原来的索引是lochead-locend
        plt_data_np = plt_data_np.transpose()#!!!需要转置


        #颜色区间0至dm
        #max=np.max(plt_data_np)
        min=0
        dm=1.5
        color=[(0,0,0)]*self.scale
        color=[color]*self.scale
        color_data=np.zeros([self.scale,self.scale])
        for i in range(self.scale):
            for j in range(self.scale):
                if plt_data_np[i][j]!=0:color_data[i][j]=(plt_data_np[i][j]-min)/dm
        rb=(10,150,0)
        rr=(0,0,255)
        for i in range(self.scale):
            color[i]=[tuple(np.add(rb,np.multiply(np.subtract(rr,rb),color_data[i][j])))for j in range(self.scale)]
        #colorbar

        font = cv2.FONT_HERSHEY_SIMPLEX  # 字体
        for i in range(self.scale):
            plt_data[i]=["%.3f"%plt_data[i][j]for j in range(self.scale)]
        text=['']*self.scale
        text=[text]*self.scale
        for i in range(self.scale):
            text[i]=[str(plt_data[i][j]) for j in range(self.scale)]

        #opencv画图模块
        for i in range(self.scale-2):
            marrayi=np.add(marray,vector2*(i+1))#每行
            for j in range(self.scale-2):
                marrayij=np.add(marrayi,vector1*(j+1))
                mpts = np.array(marrayij, np.int32)
                if (i<int(self.scale/2) and j>(int(self.scale/2)-i-2)) or (i>=int(self.scale/2) and j<=3/2*self.scale-i-4):
                    cv2.polylines(img, [mpts], True, (0, 255, 255),2)#画线
                    cv2.fillPoly(img, [mpts],color[i+1][j+1])#填充
                    #添加文字
                    r_text=[fxmid+(i+1)*d/2+(j+1)*d-5/12*d,fymid+(i+1)*l*3/4+l/8]
                    r_text=[int(r_text[0]),int(r_text[1])]
                    r_text=tuple(r_text)
                    cv2.putText(img, text[i+1][j+1], r_text, font, 0.4, (255, 255, 255), 1, cv2.LINE_AA)

        bu=df.loc['burnup', '0']
        bu=np.array(bu)

        cv2.putText(img,'RMC Power distribution at '+str(bu[bustep-1])+' MWd/kgU burnup step',(int(1/10*xmid),int(1/10*ymid)), cv2.FONT_HERSHEY_COMPLEX , int(ximg/1000), (255, 255, 255), 2, cv2.LINE_AA)
        cv2.imshow('RMC Power distribution at '+str(bu[bustep-1])+' MWd/kgU burnup step', img)

        key=cv2.waitKey(0)
        if key == 27:  # wait for ESC key to exit
            cv2.destroyAllWindows()
        elif key == ord('s'):  # wait for 's' key to save and exit
            cv2.imwrite( self.name+' burnup_step='+str(bu[bustep-1])+'.jpg', img,[int(cv2.IMWRITE_JPEG_QUALITY),100])
            cv2.destroyAllWindows()

        #plt.imshow(img, cmap='gray', interpolation='bicubic')
        #plt.xticks([-500, 500]), plt.yticks([-500, 500])  # to hide tick values on X and Y axis
        #plt.show()
        #print(L)

if __name__ == '__main__':
    dp=RMC_ANA('depletion',21)
    #dp.NPOW()
    #dp.read_NPOW()
    dp.plt_POW(7)


