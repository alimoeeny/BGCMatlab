function [ipay, left, payment] = amortize(principal, interest, term, varargin)
%[ipay, left, payment] = amortize(principal, interest, term, varargin)

monthly = 1;
payment = 0;
printdebt = 0;
bankcheat = 0;
%bankcheat mimics the way banks do it, with monthly interest = annual/12
%PNC amortization 30 years 244,000 at 3.25% = 1061.91 monthly. First month
%660.83 interest.  First year 7857.66
payed = [];
printdebt = 2;
j = 1;
while j < nargin -2
  if strncmpi(varargin{j},'bankcheat',8)
      bankcheat = 1;
  elseif strncmpi(varargin{j},'Payment',5)
        j = j + 1;
        payment = varargin{j};
  elseif strncmpi(varargin{j},'Paid',4)
      j = j+1;
      payed = varargin{j};
      payed = [0 payed];
  elseif strncmpi(varargin{j},'print',4)
      printdebt = 1;
    end
    j = j+1;
end

interest = interest/100;
if payment == 0
    intst = 1 + interest/12; %monthly interest
    if bankcheat
        intst = 1 + interest/12; %monthly interest
    else
        intst = (1+interest)^(1./12);
    end
    mterm = term * 12; %term in months
    payment = (principal * intst^(mterm))/((1-intst^(mterm))/(1-intst));
end

left(1) = principal;
isum = 0;

mi = ((1+interest)^(1/12)) -1;
if bankcheat
    mi = (interest/12);
else
    mi = ((1+interest)^(1/12)) -1;
end
fprintf('Amoritzation for $%.2f  over %d years = %d monthly payments of %.2f\n',left(1),term,term*12,payment);
if(monthly)
    fprintf('Annnual interest = %.2f%% = %.4f%% per month\n',interest*100,(intst-1) * 100);
    months = term * 12;
    if length(payed)
        months = length(payed);
    else
        months = term * 12;
        payed = ones(1,months+1).*payment;
    end
    for j = 2:months+1
        ipay(j) = mi * left(j-1);
        left(j) = left(j-1) + ipay(j) - payed(j);
        principal(j) = payed(j)-ipay(j);
    end
    payed(1) = 0;
else
    for j = 1:term+1
        ipay(j) = interest * left(j);
        ipay(j+1) = ipay(j);
        payed(j) = payment * 12;
        left(j+1) = left(j) + ipay(j) - payed(j+1);
        principal(j+1) = payed(j+1)-ipay(j);
    end
end

%calculate Yearly averages...
fprintf('Yearly totals, divided into equal monthly payments\n');
fprintf('Year\tInterest\tper month \t Principal\t per month\n');
for j = 1:term
    mid = 2+(j-1) * 12:2+j*12-1;
    isum = sum(ipay(mid));
    psum = sum(principal(mid));
    fprintf('%d    \t%.2f   \t%.2f    \t%.2f\t%.2f \n',j,isum,isum/12, psum,psum/12);
end
if printdebt
    fprintf('Full monthly amortization\n');
    last = find(left > 0);
    last = last(end);
    last = length(left);
    if printdebt == 2
        last = 2;
    end
    fprintf('Month\tInterest\tPaid      \tdebt remaining\n');
    for j = 1:last
        fprintf('%d   \t%.2f     \t%.2f    \t%.2f\n',j-1,ipay(j),payed(j),left(j));
    end
end
    end
