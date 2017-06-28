import numpy as np
from matplotlib import image as mpimg
from matplotlib import pyplot as plt
import cv2
import math
from scipy.spatial import ConvexHull


def kernel_filter(img, ordering, size):
    new_img = np.ndarray(img.shape, dtype=np.int_)    
    def apply_kernel(img, new_img, ordering, x, y, size):     
        assert size % 2 == 1
        top = max(y-(size//2), 0)
        bottom = min(y+(size//2), img.shape[0])
        left = max(x-(size//2), 0)
        right = min(x+(size//2), img.shape[1])        
        
        pixels = []        
        for i in range(top, bottom):
            for j in range(left, right):
                pixels += [img[i][j]]
        pixels = ordering(pixels)
                
        pixel_value = pixels[len(pixels)//2 + 1]
        new_img[y][x] = pixel_value

    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            apply_kernel(img, new_img, ordering, j, i, size)
    return new_img

def normal_ordering(pixels):
    ordered = sorted(pixels, key=lambda point: 0.33*point[0] + 0.33*point[1] + 0.34*point[2])
    return ordered

def marginal_ordering(pixels):    
    p_red = [p[0] for p in pixels]
    p_green = [p[1] for p in pixels]
    p_blue = [p[2] for p in pixels]
    
    if len(pixels) % 2 == 1:
        middle_red = p_red[len(p_red)//2 + 1]
        middle_green = p_green[len(p_green)//2 + 1]        
        middle_blue = p_blue[len(p_blue)//2 + 1]
    else:
        middle_red = (int(p_red[len(p_red)//2]) + p_red[len(p_red)//2 + 1]) // 2
        middle_green = (int(p_green[len(p_green)//2]) + p_green[len(p_green)//2 + 1]) // 2
        middle_blue = (int(p_blue[len(p_blue)//2]) + p_blue[len(p_blue)//2 + 1]) // 2
            
    avg_red = sum(p_red) // len(p_red)
    avg_green = sum(p_green) // len(p_green)
    avg_blue = sum(p_blue) // len(p_blue)
    
    p_red.sort()
    p_green.sort()
    p_blue.sort()
    if len(pixels) % 2 == 1:
        median_red = p_red[len(p_red)//2 + 1]
        median_green = p_green[len(p_green)//2 + 1]
        median_blue = p_blue[len(p_blue)//2 + 1]
    else:
        median_red = (int(p_red[len(p_red)//2 + 1]) + p_red[len(p_red)//2]) // 2
        median_green = (int(p_green[len(p_green)//2 + 1]) + p_green[len(p_green)//2]) // 2
        median_blue = (int(p_blue[len(p_blue)//2 + 1]) + p_blue[len(p_blue)//2]) // 2
    
    p_red.sort(key=lambda pixel: abs(int(pixel) - middle_red) + abs(int(pixel) - avg_red) + abs(int(pixel) - median_red))
    p_green.sort(key=lambda pixel: abs(int(pixel) - middle_green) + abs(int(pixel) - avg_green) + abs(int(pixel) - median_green))
    p_blue.sort(key=lambda pixel: abs(int(pixel) - middle_blue) + abs(int(pixel) - avg_blue) + abs(int(pixel) - median_blue))
    
    answer = [[p_red[i], p_green[i], p_blue[i]] for i in range(len(pixels))]
    return answer

img = mpimg.imread('../1_1.bmp')

plt.subplot(1, 3, 1)
imgplot = plt.imshow(img)
plt.title('Original Image')

plt.subplot(1, 3, 2)
filtered_img = kernel_filter(img, normal_ordering, 7)
imgplot = plt.imshow(filtered_img)
plt.title('Normal ordering')

plt.subplot(1, 3, 3)
filtered_img = kernel_filter(img, marginal_ordering, 7)
imgplot = plt.imshow(filtered_img)
plt.title('Marginal ordering')
