s=textread('strm_s.txt');
row=textread('strm_row.txt');
col=textread('strm_col.txt');
col=col';
row=row';
s=s';
omega=textread('strm_omega.txt');
nx=360;
ny=385;
keyboard
omega=reshape(omega,ny,nx);
omega=omega(:);
A=sparse(row,col,s,nx*ny,nx*ny);
strmf=reshape(A\omega,nx,ny);
strmf=strmf-strmf(25,290);
save('strm_strmf.mat','strmf')
exit
