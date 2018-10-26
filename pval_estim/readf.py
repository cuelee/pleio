def readf(afile):
    with open(afile,'r') as fin:
        n = 0;
        for line in fin:
            n = n + 1;

    with open(afile,'r') as fin:
        x = [0] * n; 
        reg = [0] * n; 
        het = [0] * n;
        i = 0;
        for line in fin:
            std = line.strip().split();
            x[i] = float(std[0]);
            reg[i] = float(std[1]);
            het[i] = float(std[2]);
            i = i + 1;
    return(x, reg, het)
