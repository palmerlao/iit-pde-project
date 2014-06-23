function [ut] = mol_1dheat_helper(t,y)
    global N D2KM KM
    coef = KM\[0;y(2:end-1);0];
    ut = 0.25.*D2KM*coef;
    ut(1)=0; ut(end)=0;
end
