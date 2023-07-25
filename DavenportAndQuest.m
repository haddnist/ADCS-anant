v1_n = [1,0,0]';
v2_n = [0,0,1]';
v1_b = [0.5,0.9,0.3]';
v1_b = v1_b/norm(v1_b);
v2_b = [0.4,0.2,0.7]';
v2_b = v2_b/norm(v2_b);
w1 = 0.1;
w2 = 1;
B = w1*v1_b*v1_n' + w2*v2_b*v2_n';
sigma = trace(B);
S=B+B';
Z = [B(2,3)-B(3,2);
     B(3,1)-B(1,3);
     B(1,2)-B(2,1)];
K = [sigma Z';
     Z (S-sigma*eye(3))];
[eigvec,eigval] = eigs(K);
maximum = max(eigval, [], 'all');
[x,y]=find(eigval==maximum);
beta_max = eigvec(:,y);
att_daven = quat2rotm(beta_max')

%quest

kappa = trace(adjoint(S));



syms x;

y = (x^2 - sigma^2 + kappa)*(x^2 - sigma^2 - (norm(Z))^2) - (x-sigma)*(Z'*S*Z+det(S))-Z'*S*S*Z;
a = w1+w2;
e = 0.0001;
N = 10000;
% Initializing step counter
step = 1;

% Finding derivate of given function
g = diff(y,x);

% Finding Functional Value
fa = eval(subs(y,x,a));

while abs(fa)> e
    fa = eval(subs(y,x,a));
    ga = eval(subs(g,x,a));
    if ga == 0
        disp('Division by zero.');
        break;
    end
    
    b = a - fa/ga;
    %fprintf('step=%d\ta=%f\tf(a)=%f\n',step,a,fa);
    a = b;
    
    if step>N
       disp('Not convergent'); 
       break;
    end
    step = step + 1;
end

fprintf('max eqigen val from quest is %f\n', a);
rho = a + sigma;
quest =[ det(rho*eye(3)-S);
        adjoint(rho*eye(3)-S)*Z
        ];
quest = quest/norm(quest);
att_quest = (quat2rotm(quest'))
