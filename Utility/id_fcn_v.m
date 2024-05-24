function cost = id_fcn_v(pars, y, u, deltaT)

    y_hat = nan(size(y));
    y_hat(1) = y(1);
   


    theta = pars(1);
    alpha(1) = pars(2);
    alpha(2) = pars(3);
    mu = pars(4);

    for i=1:(length(y)-1)
        y_hat(i+1) = y_hat(i) + (theta*(mu-y_hat(i)) + (alpha(1)*u(i,1)+alpha(2)*u(i,2)))*deltaT;
    end
        
    cost = sum((y_hat-y).^2,1);


end