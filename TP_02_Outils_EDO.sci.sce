// TP2 - EDO scalaires : Euler fourni, Heun et RK4 a completer.
// Les erreurs TODO sont volontaires : remplacer ces lignes par vos calculs.
function [t,y,nfev] = integrer_edo(f,a,b,y0,N,methode)
    if size(N,"*")<>1 then error("N doit etre scalaire"); end
    if isnan(N) | isinf(N) then error("N doit etre fini"); end
    if N<=0 | N<>floor(N) then error("N doit etre entier positif"); end
    if b<=a then error("Il faut b>a"); end
    h=(b-a)/N;
    t=a+(0:N)*h;
    y=zeros(1,N+1); y(1)=y0; nfev=0;
    for k=1:N
        tk=t(k); yk=y(k);
        select methode
       case "heun" then
    k1 = f(tk,yk);

    ypred = yk + h*k1;

    k2 = f(tk+h,ypred);

    y(k+1) = yk + h*(k1+k2)/2;

    nfev = nfev + 2;
// Heun
k1 = f(tk,yk);

predicteur = yk + h*k1;

k2 = f(tk+h,predicteur);

y(k+1) = yk + h*(k1+k2)/2;

nfev = nfev + 2;
        case "rk4" then
            // TODO R : calculer les quatre pentes et y(k+1).
            error("TODO R : completer RK4 avant de le lancer");
            nfev=nfev+4;
        else
            error("Methode inconnue : euler, heun ou rk4 attendue");
        end
    end
endfunction

function p = ordre_edo(Eh,Eh2)
    p=%nan;
    if isnan(Eh) | isnan(Eh2) | isinf(Eh) | isinf(Eh2) then return; end
    if Eh>0 & Eh2>0 then p=log(Eh/Eh2)/log(2); end
    // Un nombre fini ne garantit pas un ordre interpretable :
    // verifier aussi stabilite, regime asymptotique et arrondis.
endfunction
