function [myEBSD] = truncate_ebsd(myEBSD,xmin,xmax,ymin,ymax)
    Ebsd = myEBSD.Ebsd;
    rr = [xmin ymin (xmax-xmin) (ymax-ymin)];
    cond = inpolygon(Ebsd,rr);
    trunc_ebsd = Ebsd(cond);
    
    myEBSD.TruncEbsd = trunc_ebsd;

end