% Write out the U- and Q-matrices as source or LaTeX code to compute higher order terms.
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created June 2016, refactored December 2016
%% Write out U to inline source code for R1 and R2 with maximally mf per line
T = 7;
% s = getAsy(1,[0 1],8);
s = getAsy(0,[0 1],8);

alpha = s.alpha
mf = 2;
for k=T-1:-1:1
    ta = '        ';
    mta1 = [ta 'R1 = R1 + ('];
    lastmta1 = [ta 'R1 = R1 + 1 + ('];
    mta2 = [ta 'R2 = R2 + ('];
    if k ~= 1
        disp(['        if ( T >= ' num2str(k+1) ' )']);
        ta = '            ';
        mta1 = [ta 'R1 = R1 + ('];
        lastmta1 = [ta 'R1 = R1 + ('];
        mta2 = [ta 'R2 = R2 + ('];
    end
    mx = ceil(k/2);
    cnt = 0;
    str1 = mta1;
    str2 = mta2;
    for m = 1:mx
        str1 = [str1 num2str(s.UL(1,1,k,m),'%+.16g') '*z^(' num2str(-m) ') '];
        str2 = [str2 num2str(s.UL(1,2,k,m)*1i*4^s.alpha,'%+.16g') '*z^(' num2str(-m) ') '];
        cnt = cnt + 1;
        if (cnt == mf) || (m == mx)
            disp([str1 ')/np^' num2str(k) ';']);
            disp([str2 ')/np^' num2str(k) ';']);
            str1 = mta1;
            str2 = mta2;
            cnt = 0;
        end
    end
    mx = ceil(3*k/2);
    str1 = mta1;
    str2 = mta2;
    for m = 1:mx
        str1 = [str1 num2str(s.UR(1,1,k,m),'%+.16g') '*d^(' num2str(-m) ') '];
        str2 = [str2 num2str(s.UR(1,2,k,m)*1i*4^s.alpha,'%+.16g') '*d^(' num2str(-m) ') '];
        cnt = cnt + 1;
        if (cnt == mf) || (m == mx)
            if (k == 1) && (m == mx)
                disp([str1 ')/np^' num2str(k) ' + 1;']);
            else
                disp([str1 ')/np^' num2str(k) ';']);
            end
            disp([str2 ')/np^' num2str(k) ';']);
            str1 = mta1;
            str2 = mta2;
            cnt = 0;
        end
    end
    if k ~= 1
        disp('        end');
    end
end


%% Write out Q to inline source code for R1 and R2
T = 7;
% s = getAsy(0,[0 1],8);
s = getAsy(1,[0 1],8);

left = 1;
if left
    zp = '*z^';
    cQ = s.QL;
else
    zp = '*d^';
    cQ = s.QR;
end

alpha = s.alpha
mf = 2;
for k=T-1:-1:1
    ta = '        ';
    mta1 = [ta 'R1 = R1 + ('];
    lastmta1 = [ta 'R1 = R1 + 1 + ('];
    mta2 = [ta 'R2 = R2 + ('];
    if k ~= 1
        disp(['        if ( T >= ' num2str(k+1) ' )']);
        ta = '            ';
        mta1 = [ta 'R1 = R1 + ('];
        lastmta1 = [ta 'R1 = R1 + ('];
        mta2 = [ta 'R2 = R2 + ('];
    end
    mx = 9-k;
    cnt = 0;
    str1 = mta1;
    str2 = mta2;
    for m = mx:-1:1
        if m == 1
            str1 = [str1 num2str(cQ(1,1,k,m),'%+.16g') ' ' ];
            str2 = [str2 num2str(cQ(1,2,k,m)*1i*4^s.alpha,'%+.16g') ' ' ];
        else
            str1 = [str1 num2str(cQ(1,1,k,m),'%+.16g') zp num2str(m-1) ' ' ];
            str2 = [str2 num2str(cQ(1,2,k,m)*1i*4^s.alpha,'%+.16g') zp num2str(m-1) ' ' ];
        end
        cnt = cnt + 1;
        if (cnt == mf) || (m == 1)
            if (k == 1) && (m == 1)
                disp([str1 ')/np^' num2str(k) ' + 1;']);
            else
                disp([str1 ')/np^' num2str(k) ';']);
            end
            disp([str2 ')/np^' num2str(k) ';']);
            str1 = mta1;
            str2 = mta2;
            cnt = 0;
        end
    end
    if k ~= 1
        disp('        end');
    end
end


%% Write out QL to Latex for R1 and R2: see the asymptotic expansion in the introduction
T = 15;
s = getAsy(0,[0 1],T);
zp = '*z^';
cQ = s.QL;

alpha = s.alpha
str = '';
for k=1:T-1
    mx = 1; % 9-k; % Maximal power of z plus one
    for m = mx:-1:1
        if m == 1
            str1 = [num2str(cQ(1,1,k,m),'%.4g') ' ' ];
            str2 = [num2str(cQ(1,2,k,m)*1i*4^s.alpha,'%.4g') ' ' ];
        else
            str1 = [num2str(cQ(1,1,k,m),'%.4g') zp num2str(m-1) ' ' ];
            str2 = [num2str(cQ(1,2,k,m)*1i*4^s.alpha,'%.4g') zp num2str(m-1) ' ' ];
        end
        str = [str ' + \begin{pmatrix} ' str1 ' \\ ' str2 ' \end{pmatrix}^T \frac{1}{n^{' num2str(k) '} } '];
    end
end
disp(str)

