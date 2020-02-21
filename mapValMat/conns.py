import numpy as np

def ffcon(pre, post, input,output):
    for i in input:
        for j in output:
            pre = np.hstack([pre,i])
            post = np.hstack([post,j])
    return pre, post

def recur(pre, post, inds, selfconns=False):
    for i in inds:
        for j in inds:
            if (selfconns or ~(i==j)):
                pre = np.hstack([pre,i])
                post = np.hstack([post,j])
    return pre, post

def cullRand(pre, post):
    if(len(pre)):
        ind = np.floor(np.random.rand()*len(pre))
        pre = np.delete(pre, ind)
        post = np.delete(post, ind)
    return pre, post


def waitUntilReady(x):
    ready=False
    while(not ready):
        ready=True
        for i in range(len(x)):
            if(x[i].poll() is None):
                ready=False
