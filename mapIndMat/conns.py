import numpy as np

def ffcon(pre, post, input,output):
    for i in input:
        for j in output:
            pre = np.hstack([pre,i])
            post = np.hstack([post,j])
    return pre, post

def recur(pre, post, inds):
    for i in inds:
        for j in inds:
            if (~(i==j)):
                pre = np.hstack([pre,i])
                post = np.hstack([post,j])
    return pre, post

