# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 11:26:55 2015

@author: Aaui
"""
import numpy as np
import scipy.optimize as opt

def meanpix(im, x, y, rad):
    wid = im.shape[0]
    hei = im.shape[1]
    
    sta_x = x-rad
    sta_y = y-rad
    ste_x = 2*rad
    ste_y = 2*rad
    if sta_x < 0:
        sta_x = 0
    if sta_y < 0:
        sta_y = 0
    if sta_x + 2*rad > wid-1:
        ste_x = sta_x + 2*rad - (wid-1)
    if sta_y + 2*rad > hei-1:
        ste_y = sta_y + 2*rad - (hei-1)

    num = 0.0
    sum = np.array([0.0,0.0,0.0])
    for x in range(ste_x):
        for y in range(ste_y):
            if im[sta_x+x,sta_y+y,0]!=0.0 or im[sta_x+x,sta_y+y,1]!=0.0 or im[sta_x+x,sta_y+y,2]!=0.0:
            #if True:
                sum+=im[sta_x+x,sta_y+y]
                num+=1.0
    sum/=np.array([num,num,num])
    """
    if im[x,y,0]<sum[0] or im[x,y,1]<sum[1] or im[x,y,2]<sum[2]:
        for x in range(ste_x):
            for y in range(ste_y):
                print (im[sta_x+x,sta_y+y])
    """
    return sum

class Rbfint:
    initialized = 0;    
    #W = np.zeros((5,5))
    sig_gauss = 1.0
    
    # HARDCODED number of frequencies
    num_freq = 6.0
    
    def norm(self, x1, x2):
        #print(x1.shape, x2.shape)
        diff = (np.dot(x1-x2, self.W))
        
        #for i in range (5):
        #    diff[i]=0.0001

        """
        df = x1-x2
        if len(df[5:])!=0:
            mx = max(df[5:])
            diff = df[:6]
            diff[5] = 3.0*mx
            
        """
        
        res = np.dot(diff , diff)
        return  np.abs( res )
        
        
    def norm_diff(self, xdiff):
        #print(x1.shape, x2.shape)
        diff = (np.dot(xdiff, self.W))
        
        #for i in range (5):
        #    diff[i]=0.0001

        """
        df = x1-x2
        if len(df[5:])!=0:
            mx = max(df[5:])
            diff = df[:6]
            diff[5] = 3.0*mx
            
        """
        leng = diff.shape[0]
        res = np.zeros(leng)
        for i in range(leng):
            res[i] = - np.dot(diff[i] , diff[i])            
            #res[i] = - np.linalg.norm(diff[i])
        
        res/=(2.0*(self.sig_gauss**2))
        expres = np.exp(res)
        #res = np.dot(diff , diff)
        return expres 
        
    def g(self, dist):
        return np.exp( -dist/(2.0*(self.sig_gauss**2)) )    
    
    def __init__(self, P, Fr,Fg,Fb, sigs):
      
        self.Pts = P
        M=P.shape[0]
        #print(M)
        G=np.ones((M,M))
        self.sig_gauss=sigs[5]
        
        #print(P.shape)
        num_filt=P.shape[1]-5
        
        self.W = np.zeros((num_filt+5,num_filt+5))        
        for i in range(5):
            self.W[i,i]=sigs[i]
        for i in range(num_filt):
            j=(i)//4
            
            if j==0:
                self.W[i+5,i+5]=6.0
            elif j==1:
                self.W[i+5,i+5]=10.0 
            elif j==2:
                self.W[i+5,i+5]=10.0
            elif j==3:
                self.W[i+5,i+5]=6.0
            elif j==4:
                self.W[i+5,i+5]=10.0
            elif j==5:
                self.W[i+5,i+5]=6.0
                 
            """            
            if j==0:
                self.W[i+5,i+5]=6.0
            elif j==1:
                self.W[i+5,i+5]=8.0 
            elif j==2:
                self.W[i+5,i+5]=10.0
            elif j==3:
                self.W[i+5,i+5]=10.0
            elif j==4:
                self.W[i+5,i+5]=10.0
            elif j==5:
                self.W[i+5,i+5]=6.0
            """
            #self.W[i+5,i+5]=3.0 *(np.sqrt(np.sqrt(j+1)))
            #self.W[i,i]=5.0
        
        #print (self.W)
        
        for x in range(M):
            for y in range(x):
                dx1=P[x,:]
                dx2=P[y,:]
                cf=self.g(self.norm(dx1,dx2))
                G[x][y]=G[y][x]=cf
                #print(cf, G[x][y],G[y][x])
        
        # add STABILIZER to diagonal        
        for x in range(M):
            G[x][x]+=sigs[6]

      
        #print(G)
        
        self.A=np.empty([3,M])
        
        #A[0] is red, A[1] is green, A[2] is blue
        if True:
            self.A[0] = np.linalg.solve(G,Fr)
            self.A[1] = np.linalg.solve(G,Fg)
            self.A[2] = np.linalg.solve(G,Fb)
        else:
            self.A[0] = opt.nnls(G, Fr)[0]
            self.A[1] = opt.nnls(G, Fg)[0]
            self.A[2] = opt.nnls(G, Fb)[0]

        #thesum2 = np.dot(G, self.Ar)

        #print(self.Ar)
        #print(Fr)
        #print(thesum2)
      
        """
        for i in range(M):
            #thesum = self.DoInterp(*P[i,:])
            thesum = self.DoInterp(P[i,:])
            print( i, 'wanted %.2f %.2f %.2f  got %.2f %.2f %.2f' % (Fr[i],Fg[i],Fb[i],thesum[0],thesum[1],thesum[2]))
        """

        #altF = np.dot(G, self.A)
        #print (F, altF, end=" ")
        self.initialized = 1
    
    def DoInterp(self, newpt):
        point=np.array(newpt)
        #sum=np.zeros(3)
        #M=self.Pts.shape[0]
        
        #VECTORIZED
        diffs = self.norm_diff(self.Pts - point)
        sum2=np.dot(self.A,diffs)
        
        """
        for i in range(M):
            dx=self.Pts[i,:]
            gfunc=self.g(self.norm(dx, point))
            sum[0]+=self.A[0][i]*gfunc
            sum[1]+=self.A[1][i]*gfunc
            sum[2]+=self.A[2][i]*gfunc
            #print(i, gfunc)
        """
        
        return sum2
            
    
        