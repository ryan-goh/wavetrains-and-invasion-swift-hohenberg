% define here a stopping condition
% - v0 is old state, v1 new. Note vi=[u0; p0(iContPar)].
% - vi(1:end-1) is vector if unkowns
% - v0(end) is primary parameter
function stop = uzstop0(v1,v0)
  stop=0;
  if v1(end-1)>5; stop=1; end
end
