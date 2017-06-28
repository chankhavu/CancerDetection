import sys
import os
import cv2
import numpy as np
from matplotlib import pyplot as plt


def process_image(image_filename, working_dir):
    print('Analyzing the image %s ...' % image_filename)
    print('(1/6) Pre-processing... ', end='')
    sys.stdout.flush()
    # Read image
    img = cv2.imread(image_filename);
    gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
    
    # noise reduction and thresholding
    blured = cv2.medianBlur(gray,5)
    ret, thresh = cv2.threshold(-blured,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
     
    # Copy the thresholded image.
    im_floodfill = thresh.copy()
     
    # Mask used to flood filling.
    # Notice the size needs to be 2 pixels than the image.
    h, w = gray.shape[:2]
    mask = np.zeros((h+2, w+2), np.uint8)
     
    # Floodfill from point (0, 0)
    cv2.floodFill(im_floodfill, mask, (0,0), 255);
     
    # Invert floodfilled image
    im_floodfill_inv = cv2.bitwise_not(im_floodfill)
     
    # Combine the two images to get the foreground.
    thresh = thresh | im_floodfill_inv

    print('Completed!')
    print('(2/6) Performing morphological operations... ', end='')
    sys.stdout.flush()

    # Now, after pre-processing, we can apply morphology operations
    kernel = np.ones((9,9),np.uint8)
    opening = cv2.morphologyEx(thresh,cv2.MORPH_OPEN,kernel, iterations = 2)
    
    # sure background area
    sure_bg = cv2.dilate(opening,kernel,iterations=2)    
    border = cv2.Canny(thresh, 0, 255)
    
    # Finding sure foreground area
    dist_transform = cv2.distanceTransform(opening,cv2.DIST_L2, 5)
    ret, sure_fg = cv2.threshold(dist_transform,0.475*dist_transform.max(),255,cv2.THRESH_BINARY)
    
    # Finding unknown region
    sure_fg = np.uint8(sure_fg)
    unknown = cv2.subtract(sure_bg,sure_fg)    
    combined = dist_transform.copy()
    combined[border == 255] = dist_transform.max()
    
    print('Completed!')
    print('(3/6) Applying markers to the image... ', end='')
    sys.stdout.flush()

    # Marker labelling
    ret, markers = cv2.connectedComponents(sure_fg)
     
    # Add one to all labels so that sure background is not 0, but 1
    markers = markers+1
     
    # Now, mark the region of unknown with zero
    markers[unknown==255] = 0
    #markers[border==255] = 255
    
    img2 = cv2.bitwise_not(img, 0, 255-thresh)
    markers = cv2.watershed(img2,markers)
    
    result = cv2.medianBlur(img, 5)
    result[markers == -1] = [0,255,0]

    print('Completed!')
    print('(4/6) Getting boundaries of separated parts... ', end='')    
    sys.stdout.flush()

    # now get the boundaries of segmented cells
    left = dict([(i, 65536) for i in range(-1, markers.max()+1)])
    right = dict([(i, -1) for i in range(-1, markers.max()+1)])
    top = left.copy()
    bottom = right.copy()


    for row in range(img.shape[0]):
        for col in range(img.shape[1]):
            m = markers.item(row, col)
            if col < left[m]:
                left[m] = col
                
            if col > right[m]:
                right[m] = col
                
            if row < top[m]:
                top[m] = row
                
            if row > bottom[m]:
                bottom[m] = row

    print('Completed!')
    print('(5/6) Saving segments to separate files... ', end='')
    sys.stdout.flush()

    result2 = result.copy()
    for i in range(2, markers.max()+1):
        center = ((left[i] + right[i])//2, (top[i] + bottom[i])//2)
        
        pt1 = (left[i], top[i])    
        pt2 = (right[i], bottom[i])
        
#        l, r = max(0, center[0] - 80), min(markers.shape[1]-1, center[0] + 79)
#        t, b = max(0, center[1] - 80), min(markers.shape[0]-1, center[1] + 79)
        
        l, r = center[0]-80, center[0]+79
        t, b = center[1]-80, center[1]+79
        if l < 0:
            l, r = 0, 159
        if r > img.shape[1]:
            l, r = img.shape[1] - 160, img.shape[1]-1
        if t < 0:
            t, b = 0, 159
        if b > img.shape[0]:
            t, b = img.shape[0] - 160, img.shape[0]-1

        subimg = img[t:b+1, l:r+1].copy()
        submask = markers[t:b+1, l:r+1].copy()

        subimg[submask != i] = (255, 255, 255) 
        subimg = np.array(subimg, dtype=np.uint8)

        #subimg = np.array(img[pt1[1]:pt2[1]+1, pt1[0]:pt2[0]+1])

        filename = ('.').join(os.path.basename(image_filename).split('.')[:-1]) + '_%d.png' % i
        filename = os.path.join(working_dir, filename)
        cv2.imwrite(filename, subimg)

        cv2.rectangle(result2, pt1, pt2, (255, 0, 0))
        
    print('Completed!')
    print('(6/6) Generating diagram... ', end='')
    sys.stdout.flush()

    # Display images.
    # Show the markers    
    plt.clf()
    fig = plt.figure(figsize=(2*result2.shape[1]/100 + 1, 2*result2.shape[0]/100 + 1))
   
    subplot1 = fig.add_subplot(221)
    subplot1.imshow(thresh, 'gray')
    subplot1.set_title("Otsu Threshold and Flood Filling", size=26)
    
    subplot2 = fig.add_subplot(222)
    subplot2.imshow(combined, 'gray')
    subplot2.set_title("Distance Transformation", size=26)
    
    subplot3 = fig.add_subplot(223)
    subplot3.imshow(markers)
    subplot3.set_title("Markers (%d total)" % markers.max(), size=26)
    
    subplot4 = fig.add_subplot(224)
    subplot4.imshow(result2)
    subplot4.set_title("Added green contour to nucleuses and red boundary rectangles", size=26)
    fig.tight_layout()
 
    figname = ('.').join(os.path.basename(image_filename).split('.')[:-1]) + '_segmentation.png'
    fig.savefig(os.path.join(working_dir, figname), dpi=100)

    plt.close(fig)
    print('Completed')

    return None

if __name__ == '__main__':
    assert len(sys.argv) == 3
    process_image(sys.argv[1], sys.argv[2])
