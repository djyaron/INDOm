function res = parseAmpacFromText(arc, out)
    
    ts = textscan(arc,'%s');
    arctext = ts{1};
    
    ts = textscan(out,'%s');
    outtext = ts{1};
    
    Hf = find(ismember(arctext,'FORMATION')==1);
    [nfound,junk] = size(Hf);
    if (nfound ~= 1)
       throw(MException('parseAmpac:TextParseError',...
           ['parseAmpac: found FORMATION ',num2str(nfound),' times']));
    end
    res.Hf = str2num(arctext{Hf+2});

    % Read in cartesian coordinates
    t1 = find(ismember(outtext,'CARTESIAN')==1);
    % take last occurence of the word cartesian
    t1 = t1(size(t1,1),1);
    t1 = t1 + 7;
    done = 0;
    iatom = 1;
    while (~done)
       if (str2num(outtext{t1}) == iatom)
          element{iatom} = outtext{t1+1};
          r(1,iatom) = str2double( outtext{t1+2} );
          r(2,iatom) = str2double( outtext{t1+3} );
          r(3,iatom) = str2double( outtext{t1+4} );
          t1 = t1 + 5;
          iatom = iatom + 1;
       else
          done = 1;
          iatom = iatom -1;
       end
    end
    res.natom = iatom;
    res.r = r;
    res.element = element;

    %%
    t1 = find(ismember(outtext,'NET')==1);
    if (size(t1,1) > 0)
       % take last occurence of the word 'NET'
       t1 = t1(size(t1,1),1);
       t1 = t1 + 13;
       done = 0;
       iatom = 1;
       while (~done)
          if (str2num(outtext{t1}) == iatom)
             c_element{iatom} = outtext{t1+1};
             charge(iatom) = str2double( outtext{t1+2} );
             t1 = t1 + 4;
             iatom = iatom + 1;
          else
             done = 1;
             iatom = iatom -1;
          end
       end
       res.rho = charge;
    end

    %% Read in charges
    t1 = find(ismember(outtext,'ELECTROSTATIC')==1);
    if (size(t1,1) > 1)
       % take last occurence of the word cartesian
       t1 = t1(size(t1,1),1);
       t1 = t1 + 7;
       totalCharge = str2double( outtext{t1} );
       t1 = t1 + 7;
       done = 0;
       iatom = 1;
       while (~done)
          if (str2num(outtext{t1}) == iatom)
             c_element{iatom} = outtext{t1+1};
             charge(iatom) = str2double( outtext{t1+2} );
             t1 = t1 + 3;
             iatom = iatom + 1;
          else
             done = 1;
             iatom = iatom -1;
          end
       end
       res.esp = charge;
    end
end