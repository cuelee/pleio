import time

def cur_time(astr):
    cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime());
    print(astr + " {}".format(cur_time));

def get_time():
    atime = time.time();
    return(atime)    
