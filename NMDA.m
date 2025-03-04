function dy = NMDA(t,y,f,ft, TimeControl)
   f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t

% if (t < 0.1*TimeControl) || (t > 0.102*TimeControl)
%     Glu=0;
% end

	Rb	= 5*10^-3;     % (/uM /ms)	: binding 		
	Ru	= 12.9*10^-3;  % (/ms)	: unbinding		
	Rd	= 8.4*10^-3;   % (/ms)	: desensitization
	Rr	= 6.8*10^-3;   % (/ms)	: resensitization 
	Ro	= 46.5*10^-3;  % (/ms)	: opening
	Rc	= 73.8*10^-3;  % (/ms)	: closing

% Clements et al. 1992
	Rb	= 5*10^-3;    %(/uM /ms)	: binding 		
	Ru	= 9.5*10^-3;  %(/ms)	: unbinding		
	Rd	= 16*10^-3;   %(/ms)	: desensitization
	Rr	= 13*10^-3;   %(/ms)	: resensitization 
	Ro	= 25*10^-3;  %(/ms)	: opening
	Rc	= 59*10^-3;   %(/ms)	: closing



    rb1 = Rb * (1e3) * f;
dy = zeros(5,1);    % a column vector
NormParam = 5;
    
dy(1) =  NormParam*(Ru*y(2) -rb1*y(1));
dy(2) =  NormParam*(rb1*y(1) - (Ru+rb1)*y(2) + Ru*y(3));
dy(3) =  NormParam*(rb1*y(2) - (Ru+Rd+Ro) * y(3) + Rr*y(4) + Rc*y(5));
dy(4) =  NormParam*(Rd* y(3) - Rr * y(4));
dy(5) =  NormParam*(Ro* y(3)  - Rc * y(5));
