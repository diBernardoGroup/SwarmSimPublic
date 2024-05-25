function cost = id_cost(pars, y, u, deltaT, model)

    y_hat = nan(size(y));
    y_hat(1) = y(1);
    down_s = 10;

    theta = pars(1);
    alpha(1) = pars(2);
    alpha(2) = pars(3);
    mu = pars(4);
    
    x = y(1);
    dt = deltaT/down_s;

    for i=1:(length(y)-1)
        for j=1:down_s
            %x = x + (theta*(mu-x) + sign(x)*(alpha(1)*u(i,1)+alpha(2)*u(i,2)))*deltaT;
            [dx,~] = model(0, x, u(i,:), theta, alpha', mu);
            x = x + dx * dt;
        end
        y_hat(i+1)=x;
    end

    cost = sum((y_hat-y).^2,1) + sum((abs(y_hat)-abs(y)).^2,1);
    % cost = sum((abs(y_hat)-abs(y)).^2,1);
    % cost = sum((y_hat-y).^2,1);

end