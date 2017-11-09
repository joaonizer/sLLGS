function [N,tempo_total,ti,tf,dt]=compute_Time(gammamu0,Ms,N,tempo_total,ti)
% Computes simulation parameters from simulation configuration

if isempty(N)
    tf=tempo_total*gammamu0*Ms;
    dt=0.5;fprintf('dt = %f',dt);
    N=round((tf-ti)/dt);
    tempo_total=tf/gammamu0/Ms;
elseif isempty(tempo_total)
    dt=0.5;fprintf('dt = %f',dt);
    tf=N*dt+ti;
    tempo_total=tf/gammamu0/Ms;
elseif isempty(ti)
    erro('ti was not defined');
else
    tf=tempo_total*gammamu0*Ms;
    dt=(tf-ti)/N;
end

end