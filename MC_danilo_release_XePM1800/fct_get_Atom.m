function [Ip kappa Z Cnl l beta alphaN alphaI] = fct_get_Atom(a)

if  strcmp(a,'H1s')
        Cnl   = 2;
        Ip    = -1*13.6057./27.2113962;
        l     = 0;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = 6;
        alphaN= 4.5;
        alphaI= 0;
elseif strcmp(a,'He')
        Cnl = 2.87;
        Ip  = -1*24.59./27.2113962;
        l   = 0;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = 7;
        alphaN= 1.38;
        alphaI= 9/32;
elseif strcmp(a,'Ne')
        Cnl = 2.1;
        Ip  = -1*21.56./27.2113962;
        l   = 1; 
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = 9;
        alphaN= 2.67;
        alphaI= 1.87;
elseif strcmp(a,'Ar')
        Cnl = 2.51;
        Ip  = -1*15.76./27.2113962;
        l   = 1;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = 9.0;
        alphaN= 11.07;
        alphaI= 7.2;
elseif strcmp(a,'Kr')
        Cnl = 2.59;
        Ip  = -1*14.00./27.2113962;
        l   = 1;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = NaN;
        alphaN= 17.08;
        alphaI= 9.25;
elseif strcmp(a,'Xe')
        Cnl = 2.72;
        Ip  = -1*12.13./27.2113962;
        l   = 1;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = 9;    % 9 is a guess
        alphaN= 27.815;
        alphaI= 20.0;
elseif strcmp(a,'H2+')
        Cnl = 2.72;
        Ip  = -1*12.13./27.2113962;
        l   = 1;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = NaN;
        alphaN= NaN;
        alphaI= NaN;

elseif strcmp(a,'He1s2s-singlet')
        Cnl = 2.87;         % Cnl coefficients?
         Ip  = -1*3.9770./27.2113962;
%         Ip  = -1*(3.9770+4.7732)./(2*27.2113962);
        l   = 0;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = 7;
        alphaN= 1.38;
        alphaI= 9/32;          

elseif strcmp(a,'sodium')
        Cnl = 2.87;
         Ip  = - 4.7732./27.2113962;
%         Ip  = -1*(3.9770+4.7732)./(2*27.2113962);
        l   = 0;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = 7;
        alphaN= 1.38;
        alphaI= 9/32; 
        
        elseif strcmp(a,'He1s2s-triplet')
        Cnl = 2.87;
Ip       = -5.13908/27.212;                             %ion. pot. H  au
%         Ip  = -1*(3.9770+4.7732)./(2*27.2113962);
        l   = 0;
        kappa = sqrt(2.*abs(Ip));
        Z     = 1;
        beta  = 7;
        alphaN= 1.38;
        alphaI= 9/32; 
else
   disp('Sorry, atom is unknown.')
end
