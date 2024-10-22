% define here a stopping condition
% - v0 is old state, v1 new. Note vi=[u0; p0(iContPar)].
% - vi(1:end-1) is vector if unkowns
% - vi(end) is primary parameter
% - val is threshold for stopping
function stop = uzstop_val(v1,v0,val)
  stop=0;
  if v1(end)>val; stop=1; end
end
