% 辅助函数：平行基底
function [d_hat]=para_base(t0,t,d)
d_hat=d-dot(d,t)/(1+dot(t0,t))*(t0+t);
end