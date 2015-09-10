# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:17:15 2015

@author: Aaui
"""

import numpy as np
import scipy.misc as misc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm

from scipy import ndimage as nd
from skimage.filters import gabor_kernel

import rbfint

fullalg = 1
# (!!!) coefficients for x,y,r,g,b and gaussian sigma (and STABILIZER):
sigmas = np.array([1.0, 1.0, 2.0, 2.0, 2.0, 0.19, 0.0])
# how many runs
num_sig = 1

#init image
imgInt = mpimg.imread('img.bmp')
img = np.empty_like(imgInt, dtype=float)
#edited image
chgInt = mpimg.imread('chg.bmp')
chg = np.empty_like(chgInt, dtype=float)
#grayscale version
image = np.empty_like(chgInt[:,:,0], dtype=float)

#scale to float
imWid = img.shape[0]
imHei = img.shape[1]
for x in range(imWid):
    for y in range(imHei):
        img[x,y]=imgInt[x,y] / 255.0
        chg[x,y]=chgInt[x,y] / 255.0
        image[x,y]= (int(imgInt[x,y,0]) + int(imgInt[x,y,1]) + int(imgInt[x,y,2]))/ (3.0 * 255.0)


#gabor filters
num_freq = 4 #6
num_angl = 6
num_filt = (num_freq * num_angl)
filterval = np.zeros((imWid,imHei, num_filt)) 

if fullalg==0:
    fig2 = plt.figure(figsize=(60,27)) # 20, 50
    
angs = [2,3,5]
for j in range (num_freq):
    for i in range(num_angl):
    #0.1+0.15
    #0.05+0.075
        #ang = angs[i]
        g = gabor_kernel(0.1+0.15*j, theta= (i / 6.0) * np.pi)
        #g = gabor_kernel(1/(2.0**j), theta= (i * 0.25) * np.pi)

        filteredr = np.sqrt(nd.convolve(image, np.real(g), mode='wrap')**2 +
                   nd.convolve(image, np.imag(g), mode='wrap')**2)        
        
        filterval[:,:,(j*num_angl+i)] = filteredr  

        #if j==4:
        #    filterval[:,:,2*(j*num_angl+i)+1] = filteredr
        
        if fullalg==0:        
            b=fig2.add_subplot(num_freq,num_angl,j*num_angl+i+1)
            imgplot2 = plt.imshow(filterval[:,:,(j*num_angl+i)], cmap=cm.gray)
            #a.set_title(1/(2.0**j))
            #b.set_title(0.1+0.15*j)
            
            #b.set_title(j*num_angl+i)
            b.set_title("% Frequency = %.2f, Direction = %d pi / 6" % (j*num_angl+i, (0.1+0.15*j), i))

#print (filterval)       
#nrmfilt = np.linalg.norm(filterval)
#filterval/=nrmfilt
#print(np.linalg.norm(filterval))

#figure with resulting images        
fig = plt.figure(figsize=(20,20)) # 20, 50
a=fig.add_subplot(2+num_sig,2,1)
imgplot = plt.imshow(chg)
a.set_title('Before')            
        
#compare images, specify edits
whiteF = []
blackF = []
resWF = []
resBF = []

#how close to the resulting color
divd=1

#684x912
for x in range(imWid):
    for y in range(imHei):
        if chg[x,y,0]!=img[x,y,0] and chg[x,y,1]!=img[x,y,1] and chg[x,y,2]!=img[x,y,2]: #modified pixels only      
            newpt = np.zeros(5+num_filt+3) # xyrgb, filters, constraints
            newpt[0:5] = np.array( [x/float(imWid),y/float(imHei),img[x,y,0],img[x,y,1],img[x,y,2]] ) # increase size
            newpt[5:(num_filt+5)] = filterval[x,y]  
            
            if chg[x,y,0]==0 and chg[x,y,1]==0 and chg[x,y,2]==0:  # marked as black
                #newpt.extend([0, 0, 0])
                newpt[(num_filt+5):]=np.array( [0.,0.,0.] )
                blackF.append(newpt) 
            else: # marked as white
                newpt[(num_filt+5):]=np.array( [(chg[x,y,0]-img[x,y,0])/divd, (chg[x,y,1]-img[x,y,1])/divd, (chg[x,y,2]-img[x,y,2])/divd]  )
                whiteF.append(newpt) 

print("Total points: ", len(whiteF) + len(blackF))

"""
diff = np.zeros(num_filt)
diffs_med = np.zeros((len(whiteF),5+num_filt+3))
for i in range(len(whiteF)):
    sum=np.zeros(5+num_filt+3)
    for j in range(len(blackF)):
        sum+=abs(whiteF[i]-blackF[j])
    sum/=len(blackF)
    diffs_med[i]=sum

for j in range(num_filt):        
    max=0.0
    for i in range(len(whiteF)):
        if diffs_med[i,j+5]>max:
            max=diffs_med[i,j+5]
    diff[j]=max
        
print(diff)
"""

#filter out the majority of the points
filt = (len(whiteF) + len(blackF))//150 # how many to filter - 50
white = []
white.extend(whiteF[1::(2*filt)])
black = []
#print(white)
black.extend(blackF[1::(filt)])
white.extend(black)

#pts = all the points for rbf
pts = np.array(white)
print("Used point: ", len(pts))
#print(pts)

endCoord = 5+num_filt 

if True:
    a=fig.add_subplot(2+num_sig,2,2)
    img2 = np.copy(img)
    img2gab = np.copy(img)
    img2diff = np.copy(img)
    img2diffinit = np.copy(img)

ims = [img2]
imsGab = [img2gab]
imsDif = [img2diff]
imsDifInit = [img2diffinit]


for n in range(num_sig):
    print (n)
    
    rbfTst = rbfint.Rbfint(pts[:,:5], pts[:,endCoord],pts[:,endCoord+1],pts[:,endCoord+2], sigmas) # 5,6,7  
    rbfGabor = rbfint.Rbfint(pts[:,:endCoord], pts[:,endCoord],pts[:,endCoord+1],pts[:,endCoord+2], sigmas) # 5,6,7

    imgN = ims[n]
    imgNgab = imsGab[n]
    imgNdiff = imsDif[n]
    imgNdiffInit = imsDifInit[n]
    
    thresh = 0.2    
    for x in range(imWid):
        for y in range(imHei):
            newpt = [x/float(imWid),y/float(imHei),img[x,y,0],img[x,y,1],img[x,y,2]]
            
            if True:
                pixel=np.array(rbfTst.DoInterp(newpt))        
                imgN[x,y]+=pixel 
                imgNdiffInit[x,y] = np.array([ 0.5,0.5,0.5 ]) + pixel
            
            for i in range(num_filt): 
                newpt.extend([filterval[x,y,i]])
                
            if True:
                pixelgab = np.array(rbfGabor.DoInterp(newpt)) 
                
                """
                if (np.abs(pixelgab[0])>thresh or np.abs(pixelgab[1])>thresh or np.abs(pixelgab[2])>thresh):
                    for k in range(3):
                        if (np.abs(pixelgab[k])>thresh):
                            newpixg = (np.sign(pixelgab[k]))*(np.sqrt(np.abs(pixelgab[k])-thresh))
                            pixelgab[k]=newpixg
                else:
                    pixelgab=np.array([0.0,0.0,0.0])
                """                    
                        
                #if (np.abs(pixelgab[0])>thresh or np.abs(pixelgab[1])>thresh or np.abs(pixelgab[2])>thresh):
                if True:
                    imgNgab[x,y]+=pixelgab
                    imgNdiff[x,y] = np.array([ 0.5,0.5,0.5 ]) + pixelgab
    

    """
    imgNblur = np.empty_like(imgNdiff, dtype=float)
    for x in range(imWid):
        for y in range(imHei):
            imgNblur[x,y]=rbfint.meanpix(imgNdiff, x, y, 1)
    imgNgab +=  imgNblur
    """
    
    imgN[np.where(imgN>1.)] = 1.
    imgN[np.where(imgN<0.)] = 0.       
        
    imgNgab[np.where(imgNgab>1.)] = 1.
    imgNgab[np.where(imgNgab<0.)] = 0.          
                 
    if n==0: 
        imgplot = plt.imshow(imgNdiff) 
        a.set_title('With Gabor')
   
    a=fig.add_subplot(2+num_sig,2,2*(n+1)+1)
    imgplot = plt.imshow(imgN) 
    misc.imsave('res1.png',img2)   
    misc.imsave('dif1.png',imgNdiffInit)  

    a.set_title('sig = %f' % sigmas[5])
    a=fig.add_subplot(2+num_sig,2,2*(n+1)+2)
    imgplot = plt.imshow(imgNgab)
    misc.imsave('res2.png',img2gab)    
    misc.imsave('dif2.png',imgNdiff)  
    
    ims.append(np.copy(img))
    imsGab.append(np.copy(img))
    imsDif.append(np.copy(img))
    
    #increase the STABILIZER
    sigmas[6]+=1.0
    











