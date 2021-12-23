function [Design] = DCE_DESIGN(Design)

    fprintf('Komentarz: Najpierw drukowane sa informacje o designie\n')
    fprintf('Komentarz: W rzeczywistosci tworzone sa potrzebne elementy \n')
    fprintf('Komentarz: Wyliczane ilosci poziomow, sprawdzane dane \n')
    % type, b1... -> check
    
    if ~isfield(Design,"number_of_alternatives")
        cprintf(rgb('Black'),'Missing element: Design.number_of_alternatives\n')
        return
    end
    if ~isfield(Design,"number_of_attributes")
        cprintf(rgb('Black'),'Missing element: Design.number_of_attributes\n')
        return
    end    
    if Design.number_of_alternatives <= 0
        cprintf(rgb('DarkOrange'),'WARNING:  Number of alternatives cannot be smaller than 0\n')
        return
    end
    if Design.number_of_attributes <= 0
        cprintf(rgb('DarkOrange'),'WARNING:  Number of attributes cannot be smaller than 0\n')
        return
    end    
    
    nattr_aux = Design.number_of_attributes;
    nalt_aux = Design.number_of_alternatives;
    
    fprintf('DESIGN \nDesign for Discrete Choice Experiment \n\n');
    fprintf('Number of alternatives : %d \n',Design.number_of_alternatives);
    fprintf('Number of attributes   : %d \n',Design.number_of_attributes);
    
    
    % is (stated) number of attributes and alternatives correct?
    for i=1:nalt_aux
        nalt_aux_for = i;
        if ~isfield(Design,strcat('Alt_name_',num2str(i)))
            fprintf('\nMissing alternative name: Design.%s\n',strcat('Alt_name_',num2str(i)));
            return
        end
        for i=1:nattr_aux
            if ~isfield(Design,strcat('a',num2str(i),'_name'))
                fprintf('\nMissing attribute name: Design.%s \n',strcat('a',num2str(i),'_name'));
                return
            end
            
            if ~isfield(Design,strcat('a',num2str(i),'_',num2str(nalt_aux_for)))
                if ~isfield(Design,strcat('a',num2str(i)))
                    fprintf('\nMissing element): a%d_%s', i, num2str(nalt_aux_for))
                    fprintf('\nor a%d (if you want to keep attribute fixed across alternatives)', i)
                    fprintf('\nEventually set value to [missing] for alternatives that does not include specific attributes\n')
                    return
                else
                    Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))) = Design.(strcat('a',num2str(i)));
                end
            end
            
            if ~isfield(Design,strcat('a',num2str(i),'_',num2str(nalt_aux_for)))
                fprintf('\nMissing element: Design.%s \nRemember to create missing fields (=[missing]) for alternatives that does not include specific attributes\n',strcat('a',num2str(i),'_',num2str(nalt_aux_for)));
                return
            end

            if ~isfield(Design,strcat('b',num2str(i),'_',num2str(nalt_aux_for)))
                if ~isfield(Design,strcat('b',num2str(i)))
                    fprintf('\nMissing element): b%d_%s', i, num2str(nalt_aux_for))
                    fprintf('\nor b%d (if you want to keep parameter fixed across alternatives)', i)
                    fprintf('\nEventually set value to [missing] for alternatives that does not include specific attributes\n')
                    return
                else
                    Design.(strcat('b',num2str(i),'_',num2str(nalt_aux_for))) = Design.(strcat('b',num2str(i)));
                end
            end
            
            if Design.(strcat('a',num2str(i),'_type'))==["D"]
                if size(Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))),2)~=size(Design.(strcat('b',num2str(i),'_',num2str(nalt_aux_for))),2)
                    fprintf('\nDiscrepancy between a and b: %d_%s', i, num2str(nalt_aux_for))
                    return
                end
            end
        end
    end
    
    ff_ncs = 1;
    for i=1:nalt_aux
        nalt_aux_for = i;
        for i=1:nattr_aux
            ff_ncs=ff_ncs*double(size(Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))),2));
        end
    end
    
    if Design.FF == 1
        Design.full_factorial = [missing];
        Design.ff_number_of_cs = ff_ncs;
        fprintf('Number of choice sets in full factorial   : %d \n \n',Design.ff_number_of_cs);
        if ff_ncs>1000000 && Design.ff_number_of_cs<1000000000
            cprintf(rgb('LightBlue'),'Remember:  Bigger amount of alternatives and attributes might increase computation time\n')
        elseif Design.ff_number_of_cs>=1000000000
            cprintf(rgb('LightBlue'),'Remember:  Bigger amount of alternatives and attributes might increase computation time\n Be patient ... it may take a while \n')
        end
    end
    
    cprintf(rgb('Black'),'\n--------------------------------------------------------------------------------\n')
    
    
    attr_type = ["Aux"];
    alt_number = ["Aux"];
    parameters = ["Aux"];
    
    vector_count = 1;
    
    for i=1:nalt_aux
        nalt_aux_for = i;

        fprintf('\n ---------- Alternative %d: %s ---------- \n',i,Design.(strcat('Alt_name_',num2str(i))));
        for i=1:nattr_aux
            nattr_aux_for = i;
            
            fprintf('Attribute %d \n',i);
            fprintf('     -> Attribute name: %s \n',Design.(strcat('a',num2str(i),'_name')));
            fprintf('     -> Attribute type: %s \n',Design.(strcat('a',num2str(i),'_type')));

            if ismissing(Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))))
                fprintf(2,'     <<<This alternative does not include this attribute>>> \n');

            else
                attr_type(end+1) = Design.(strcat('a',num2str(i),'_type'));
                alt_number(end+1) = nalt_aux_for;
                fprintf('     -> Number of levels: %d \n',double(size(Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))),2)));
              
                if Design.FF == 1
                    Design.full_factorial = combvec(Design.full_factorial,vector_count:vector_count-1+double(size(Design.(strcat('a',num2str(nattr_aux_for),'_',num2str(nalt_aux_for))),2)));
                end
               
                for i=1:double(size(Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))),2))

                    fprintf('          # Level %d: %s \n',i,Design.(strcat('a',num2str(nattr_aux_for),'_',num2str(nalt_aux_for)))(i));

                    vector_name(vector_count)=string(Design.(strcat('a',num2str(nattr_aux_for),'_',num2str(nalt_aux_for)))(i));
                    vector_count = vector_count+1;
                    if Design.(strcat('a',num2str(nattr_aux_for),'_type'))=="D"
                        parameters(end+1) = Design.(strcat('b',num2str(nattr_aux_for),'_',num2str(nalt_aux_for)))(i);
                    else
                            parameters(end+1) = Design.(strcat('b',num2str(nattr_aux_for),'_',num2str(nalt_aux_for)));
                    end
                end
            end
        end  
    end

    size_vector_name = size(vector_name,2);
    if Design.FF == 1
        Design.full_factorial = strcat(string((Design.full_factorial)'),'*');
        Design.full_factorial(:,1)=[];
    end
    
    attr_type(:,1)=[];
    alt_number(:,1)=[];
    parameters(:,1)=[];
    alt_number = double(alt_number);
    
    
    cprintf(rgb('Black'),'\n--------------------------------------------------------------------------------\n')
    
    vector_count_random = 1;
    
    % Random choice sets
    if isfield(Design,"RandomCS")
        fprintf('Number of random choice sets (random fractional factorial): %d\n', Design.RandomCS')
        Design.fractional_factorial_random = [];
        if Design.FF == 1
            randomff = randperm(size(Design.full_factorial,1),Design.RandomCS);       
            for i=1:Design.RandomCS
                Design.fractional_factorial_random = [Design.fractional_factorial_random;Design.full_factorial(randomff(i),:)];
            end
        else
            for i=1:nalt_aux
                nalt_aux_for = i;
                for i=1:nattr_aux
                    if ~ismissing(Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))))
                        Design.fractional_factorial_random(:,end+1) = randsample(vector_count_random:vector_count_random-1+double(size(Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))),2)),Design.RandomCS,1);
                        vector_count_random = vector_count_random + double(size(Design.(strcat('a',num2str(i),'_',num2str(nalt_aux_for))),2));
                    end
                end
            end
            Design.fractional_factorial_random = strcat(string((Design.fractional_factorial_random)),'*');
        end
    end      
   
   if Design.EF_mnl==1
       precision = Design.EF_prec;
       simulation = Design.EF_sim;
       numbercs = Design.EF_cs;
       min_cond = Design.EF_error;
       Design.FisherIN = [];
       Design.UB = [0];
       UB_aux = Design.UB;
       Design.Derror = [999*999e999];
       Derror_aux = Design.Derror;
       Design.Aerror = [999*999e999];
       Aerror_aux = Design.Aerror;
       Design.Sestimate = [999*999e999];
       Sestimate_aux = Design.Sestimate;

       % Efficient Design
       for h=1:simulation
            fprintf('Komentarz: Czesc, oto prezentacja jak wyglada aktualnie kod wyszukujacy najlepszy efficient Design\n')
            fprintf('Komentarz: Kod nie jest dokonczony. Jeszcze planuje wiele usprawnien i poprawek\n')
            fprintf('Komentarz: Na ta chwile nie wiem dlaczego czasem w macierzy AVC pojawiaja sie ujemne wartosci\n')
            fprintf('Komentarz: Mam nadzieje, ze ta prezentacja przyblizy Ci to na czym stoje i moze wylapiesz miejsce gdzie lezy blad\n')
            % losowanie
            fractional_factorial_random = [];
            vector_count_random = 1;

            fprintf('\nZaczynajmy :)\n')
            fprintf('Komentarz: Najpierw losowowany jest design (random fractional factorial)\n')
            fprintf('Komentarz: Jesli wybierze sie opcje Design.EF=1 to najpierw tworzony jest full factorial\n')
            fprintf('Komentarz: Oczywiscie warto pamietac, ze full factorial moze byc zbyt duzy, jesli jest duzo atrybutow, poziomow lub alternatyw\n')
            fprintf('Komentarz: Jesli nie ma utworzonego full factoriala to design jest losowany na podstawie losowo wybranych poziomow\n')
            if Design.FF == 1
                randomff = randperm(size(Design.full_factorial,1),numbercs);       
                r=1;
                for r=1:numbercs
                    fractional_factorial_random = [fractional_factorial_random;Design.full_factorial(randomff(r),:)];
                end
            else
                q=1;
                for q=1:nalt_aux
                    nalt_aux_for = q;
                    a=1;
                    for a=1:nattr_aux
                        if ~ismissing(Design.(strcat('a',num2str(a),'_',num2str(nalt_aux_for))))
                            fractional_factorial_random(:,end+1) = randsample(vector_count_random:vector_count_random-1+double(size(Design.(strcat('a',num2str(a),'_',num2str(nalt_aux_for))),2)),numbercs,1);
                            vector_count_random = vector_count_random + double(size(Design.(strcat('a',num2str(a),'_',num2str(nalt_aux_for))),2));
                        end
                    end
                end
                fractional_factorial_random = strcat(string((fractional_factorial_random)),'*');
            end

            fprintf('Komentarz: Tak wyglada wylosowany design\n')
            fractional_factorial_random
            fprintf('Komentarz: Te gwiazdki sa technicznie potrzebne, by odroznic wartosci liczbowe od atrybutow\n')
            fprintf('Komentarz: W pozniejszej fazie kodu na podstawie tych informacji z gwiazdkami odnajdowane sa wartosci priorow \nlub informacje czy dana zmienna jest ciagla, badz dyskretna\n')
            attributes_without_star = double(erase(fractional_factorial_random,"*"));
            
            aux=1;
            decoded_parameters = [];
            alt_number_decoded = [];
            attr_type_decoded = ["Aux"];
            min_aux=[];
            max_aux=[];

            fprintf('Komentarz: Zmienne dyskretne poczatkowo sa w jednej kolumnie\n')
            fprintf('Komentarz: Potrzeba je jednak rozkodowac na poszczegolne zmienne binarne\n')
            iter=1;
            for iter=1:nalt_aux
                nalt_aux_for = iter;
                uiter=1;
                for uiter=1:nattr_aux
                    if Design.(strcat('a',num2str(uiter),'_type'))=="D"
                        if ~ismissing(Design.(strcat('a',num2str(uiter),'_',num2str(nalt_aux_for))))
                            min_aux = min(attributes_without_star(:,aux));
                            max_aux = max(attributes_without_star(:,aux));
                            diter=min_aux;
                            for diter=min_aux:max_aux
                                cond = (attributes_without_star(:,aux) == diter);
                                decoded_parameters(cond,end+1)=attributes_without_star(cond,aux);
                                alt_number_decoded(end+1)=nalt_aux_for;
                                attr_type_decoded(end+1)=Design.(strcat('a',num2str(uiter),'_type'));
                            end
                            aux = aux + 1;
                        end
                    else    
                        if ~ismissing(Design.(strcat('a',num2str(uiter),'_',num2str(nalt_aux_for))))
                            decoded_parameters(:,end+1)=attributes_without_star(:,aux);
                            alt_number_decoded(end+1)=nalt_aux_for;
                            attr_type_decoded(end+1)=Design.(strcat('a',num2str(uiter),'_type'));
                            aux = aux + 1;
                        end                 
                    end
                end
            end

            fprintf('Komentarz: Tak wyglada rozkodowany design\n')
            decoded_parameters
            fprintf('Komentarz: Przy okazji stworzony zostal pomocniczy wektor informujacy o numerze alternatywy danego atrybutu\n')
            alt_number_decoded
            fprintf('Komentarz: Oraz rodzaju atrybutu\n')
            
            decoded_parameters = strcat(string(decoded_parameters),'*');
            attr_type_decoded(:,1)=[]


            parametersFI = decoded_parameters;
            attributesFI = decoded_parameters;

            kiter=1;
            for kiter=1:size_vector_name
                parametersFI(parametersFI==strcat(num2str(kiter),'*'))=parameters(kiter);
                attributesFI(attributesFI==strcat(num2str(kiter),'*'))=vector_name(kiter);
            end    

            fprintf('Komentarz: Teraz odnajdowane sa odpowienie nazwy atrybutow oraz parametry (do tego byly potrzebne gwiazdki)\n')
            fprintf('Komentarz: Parametry:\n')
            parametersFI
            fprintf('Komentarz: Oraz atrybuty\n')
            attributesFI
            
            parametersFI(attr_type_decoded=="D" & decoded_parameters=="0*")=0;
            parametersFI = double(parametersFI);


            attr_type_decoded = repmat(attr_type_decoded,size(fractional_factorial_random,1),1);
            attributesFI(attr_type_decoded=="D" & decoded_parameters~="0*")=1;
            attributesFI(attr_type_decoded=="D" & decoded_parameters=="0*")=0;

            attributesFI=double(attributesFI);

            attributes_final = [0];
            attributes_final = repmat(attributes_final,size(attributesFI,1),1);
            parameters_final = [0];
            parameters_final = repmat(parameters_final,size(attributesFI,1),1);
            alt_number_decoded2 = [];

            fprintf('Komentarz: Po malych transofmacjach:\n')
            fprintf('Komentarz: Atrybuty\n')
            attributesFI            
            fprintf('Komentarz: Parametry\n')
            parametersFI

            
            fprintf('Komentarz: brakuje jednak stalych\n')
            % adding constants
            jiter=1;
            for jiter=1:nalt_aux
                value_aux = jiter;
                if ~ismissing(Design.(strcat('b0_',num2str(jiter))))
                    parameters_final(:,end+1)=Design.(strcat('b0_',num2str(jiter)));
                    attributes_final(:,end+1)=1;
                    alt_number_decoded2(:,end+1)=value_aux;

                end
                fiter=1;
                for fiter=1:size(attributesFI,2)
                    if alt_number_decoded(fiter)==value_aux
                        attributes_final(:,end+1)=attributesFI(:,fiter);
                        parameters_final(:,end+1)=parametersFI(:,fiter);
                        alt_number_decoded2(:,end+1)=alt_number_decoded(:,fiter);
                    end
                end
            end
            attributes_final(:,1)=[];
            parameters_final(:,1)=[];
            parameters_nan = isnan(sum(parameters_final,1));

            fprintf('Komentarz: Po dodaniu stalych\n')
            fprintf('Komentarz: Atrybuty:\n')
            attributes_final
            fprintf('Komentarz: Parametry:\n')
            parameters_final
            fprintf('Komentarz: Oraz poprawione numery alternatyw:\n')
            alt_number_decoded2

            % deleting columns with NaN (base levels)
            oiter=1;
            for oiter=1:size(parameters_final,2)
                if parameters_nan(oiter)==1
                    attributes_final(:,oiter)=[NaN];
                    parameters_final(:,oiter)=[NaN];
                    alt_number_decoded2(:,oiter)=[NaN];

                end
            end


            attributes_final = attributes_final(:,~all(isnan(attributes_final)));
            parameters_final = parameters_final(:,~all(isnan(parameters_final)));
            alt_number_decoded2 = alt_number_decoded2(~isnan(alt_number_decoded2));

            fprintf('Komentarz: Przy okazji mozemy sie pozbyc kolumn dot. poziomow bazowych (dla zmiennych dyskretnych)\n')
            fprintf('Komentarz: Atrybuty\n')
            attributes_final
            fprintf('Komentarz: Parametry\n')
            parameters_final
            fprintf('Komentarz: Numery alternatyw \n')
            alt_number_decoded2
            
            
            fprintf('Komentarz: Majac przygotowane w odpowiednim formacie dane wejsciowe mozemy zabrac sie do wyliczania Fischer Information matrix\n')
            fprintf('Komentarz: Macierz Fishera to macierz o wymiarach A na A, gdzie A to wektor uzytych parametrow\n')
            fprintf('Komentarz: Macierz Fishera to obliczenia pochodnych dla funkcji LL wzgledem parametrow\n')
            fprintf('Komentarz: Ponizej przyklad wygladu macierzy Fishera dla 4 parametrow\n')
            ["b1b1","b1b2","b1b3","b1b4";"b2b1","b2b2","b2b3","b2b4";"b3b1","b3b2","b3b3","b3b4";"b4b1","b4b2","b4b3","b4b4"]
            fprintf('Komentarz: Informuje to o tym, ktore pochodne trzeba wyliczyc\n')
            fprintf('Komentarz: Np. b3b4 -> pochodna LL po b3 a nastepnie po b4\n')
            fprintf('Komentarz: Mala uwaga: jest to macierz symetryczna b3b1=b1b3\n')
            fprintf('\nKomentarz: Jak przyblizyc pochodna drugiego stopnia?\n')
            fprintf('Komentarz: Skoro znamy priory to ruszmy odpowiednie parametry o mala ilosc i sprawdzmy jak zmienia sie dla nich funkca LL\n')
            fprintf('Komentarz: Pochodna niesie bowiem informacje o nachyleniu danej funkcji (http://matlablog.ont.com.pl/rozniczkowanie-numeryczne-czyli-nie-taka-straszna-pochodna/)\n')
            fprintf('Komentarz: Tu wzor (dla przyblizenia) pochodnej drugiego stopnia dla dwoch parametrow:\n')
            fprintf('Komentarz: Multivariate finite differences (na samym dole) (https://en.wikipedia.org/wiki/Finite_difference)\n')
            fprintf('Komentarz: Potrzeba jedynie wyliczyc LL dla funkcji w ktorej ruszamy parametrami i wstawic do odpowiedniego wzoru\n')
            fprintf('Komentarz: Na przekatnej mozemy uzyc wzoru fxx(x,y)~=(f(x+h,y)-2*f(x,y)+f(x-h,y))/(h^2) \n')
            fprintf('Komentarz: A poza przekatna: fxy(x,y)~=(f(x+h,y+k)-f(x+h,y-k)-f(x-h,y+k)+f(x-h,y-k))/(4*h*k) \n')
            fprintf('Komentarz: Kod dla kazdego przypadku wyglada tak samo\n')
            fprintf('Komentarz: np. LL_plusplus odpowiada za wyliczenie f(x+h,y+k) i dla tego przypadku oba parametry sa powiekszane\n')
            fprintf('Komentarz: LL_plusminus odpowiada za f(x+h,y-k); LL_minusplus za f(x-h,y+k); LL_minusminus za f(x-h,y-k) \n')
            fprintf('Komentarz: Przygotowane wczesniej parametry wystarczy odpowiednio zmienic, a nastepnie wyliczyc dla nich LL\n')
            fprintf('\nKomentarz: Ponizej przyklad dla LL bez zmian parametrow - potrzebne do wzoru dla przekatnej\n')
            
            % Fisher Information Matrix

            fprintf('Komentarz: Najpierw usuwamy kolumny z brakiem danych (glownie dotyczy do zmiennych dyskretnych)\n')
            fprintf('Komentarz: Jesli wylosowany zestaw choice taskow bez jednego z poziomow dyskretnych to nie powinien on byc wyliczany\n')
            fprintf('Komentarz: Bety dla tego parametru nie ma (bo ten poziom nie wystepuje), wiec nie mozna dla tego przypadku obliczyc macierzy FI\n')
            % deleting columns with zeros (absent attribute)
            aux_absent = find(all(attributes_final==0));        
            attributes_final(:,aux_absent)=[];
            parameters_final(:,aux_absent)=[];
            alt_number_decoded2(:,aux_absent)=[];

            fprintf('Komentarz: Dla przypomnienia\n')
            fprintf('Komentarz: Wylosowany design: (tu nie uzywany -> gdyz mamy przygotowane odpowiednie atrybuty)\n')
            fractional_factorial_random
            fprintf('Komentarz: Przygotowane atrybuty:\n')
            attributes_final
            fprintf('Komentarz: Przygotowane parametry:\n')
            parameters_final

            % LL - without change

            fprintf('Komentarz: Wyliczamy V (uzytecznosc  dla kazdego parametru)\n')
            fprintf('Komentarz: V = attributes_final.*parameters_final \n')
            V = attributes_final.*parameters_final

            LL=[];
            LL_max=[];
            iy=1;
            for iy=1:size(parameters_final,1)
                aux_LL = iy;
                for dp=1:nalt_aux
                    LL(aux_LL,dp)=sum(V(aux_LL,alt_number_decoded2==dp));
                    LL_max(aux_LL,1)=max(LL(aux_LL,:));
                end
            end
            
            fprintf('Komentarz: Oczywiscie interesuje nas LL dla oddzielnych alternatyw\n')
            fprintf('Komentarz: Chcemy zatem zsumowac czastkowe uzytecznosci dla kazdej z alternatyw\n')
            LL
            fprintf('Komentarz: Przy okazji znajdzmy alternatywe z najwieksza wartoscia LL -> ta alternatywa zostanie "wybrana" \n')
            LL_max
            
    %         LL
    %         LL_max

            fprintf('Komentarz: Nalezy odnalezc teraz wybrana alternatywe (potrzebne do wyliczenia prawdopodobienstwa wyboru alternatyw)\n')
            % finding max LL
            ij=1;
            choice_aux=[];
            Choice=zeros(size(parameters_final,1),nalt_aux);
            for ij=1:size(parameters_final,1)
                choice_aux(ij,1)=find(LL(ij,:)==LL_max(ij,1),1,'first');
                Choice(ij,choice_aux(ij,1))=1;
            end

            fprintf('Komentarz: Funkcja odpowiednio znajduje polozenie LL_max w LL \n')
            Choice
            fprintf('Komentarz: Nie mozna po prostu znalezc LL_max, gdyz istnieje przypadek dla rownomiernie rozlozonego prawdopodobienstwa\n')
            fprintf('Komentarz: Wybrana alternatywa moze byc tylko 1\n')
            fprintf('Komentarz: Na szczescie nie interesuje nas dokladnie ktora jest wybrana (jesli sa dwie o tym samym LL)\n')
            fprintf('Komentarz: Prawdopodobienstwo wyboru takich alternatyw jest identyczne, wiec dla nas obojetnie ktora wejdzie do rownania\n')
            fprintf('Komentarz: Prob = exp(V_wybrana)/suma{exp(V_i)}\n')

            fprintf('Komentarz: Najpierw bierzemy exp(LL)\n')
            LL=exp(LL)
            fprintf('Komentarz: Wyliczamy wyrazenie na gorze LL_up=LL.*Choice\n')
            LL_up = LL.*Choice;

            fprintf('Komentarz: Wyliczamy prawdopodobienstwa LLprob(i,1) = sum(LL_up(i,:))/sum(LL(i,:))')
            is=1;
            LLprob=[];
            for is= 1:size(parameters_final,1)
                LLprob(is,1) = sum(LL_up(is,:))/sum(LL(is,:));
            end

            LLprob
            fprintf('Komentarz: Zlogarytmujmy')
            LLprob = log(LLprob)
            fprintf('Komentarz: I liczymy sume (calkowite LL)')
            LL_value = sum(LLprob);

            
            fprintf('Komentarz: Analogicznie wyliczone sa potrzebne LL_plusplus etc. \n')
            fprintf('Komentarz: (odpowiednio dodajemy/odejmujemy zadeklarowana wartosc dla utworzonej macierzy z parametrami)\n')
            
            FisherIN = [];

            % parameterplus (only for =~ 0)
            % diagonal
            yf=1;
            for yf=1:size(parameters_final,2)
                parameters_plus = parameters_final;
                parameters_plus(:,yf) = parameters_plus(:,yf)+precision;

                parameters_plusplus=parameters_final;
                parameters_plusplus(:,yf) = parameters_plusplus(:,yf)+(2*precision);

                V_plus = attributes_final.*parameters_plus;
                V_plusplus = attributes_final.*parameters_plusplus;

                tl=1;
                LL_plus = [];
                LL_plus_max = [];
                LL_plusplus = [];
                LL_plusplus_max = [];
                for tl=1:size(parameters_final,1)
                    aux_LL = tl;
                    sk=1;
                    for sk=1:nalt_aux
                        LL_plus(aux_LL,sk)=sum(V_plus(aux_LL,alt_number_decoded2==sk));
                        LL_plus_max(aux_LL,1)=max(LL_plus(aux_LL,:));

                        LL_plusplus(aux_LL,sk)=sum(V_plusplus(aux_LL,alt_number_decoded2==sk));
                        LL_plusplus_max(aux_LL,1)=max(LL_plusplus(aux_LL,:));
                    end
                end

                % finding max LL
                gl=1;    
                choice_plus_aux = [];
                Choice_plus = zeros(size(parameters_final,1),nalt_aux);
                choice_plusplus_aux = [];
                Choice_plusplus = zeros(size(parameters_final,1),nalt_aux);
                for gl=1:size(parameters_final,1)
                    choice_plus_aux(gl,1)=find(LL_plus(gl,:)==LL_plus_max(gl,1),1,'first');
                    Choice_plus(gl,choice_plus_aux(gl,1))=1;

                    choice_plusplus_aux(gl,1)=find(LL_plusplus(gl,:)==LL_plusplus_max(gl,1),1,'first');
                    Choice_plusplus(gl,choice_plusplus_aux(gl,1))=1;
                end

                LL_plus=exp(LL_plus);
                LL_plus_up = LL_plus.*Choice_plus;

                LL_plusplus=exp(LL_plusplus);
                LL_plusplus_up = LL_plusplus.*Choice_plusplus;

                cf=1;
                LLprob_plus = [];
                LLprob_plusplus = [];
                for cf= 1:size(parameters_final,1)
                    LLprob_plus(cf,1) = sum(LL_plus_up(cf,:))/sum(LL_plus(cf,:));
                    LLprob_plusplus(cf,1) = sum(LL_plusplus_up(cf,:))/sum(LL_plusplus(cf,:));
                end


                LLprob_plus = log(LLprob_plus);
                LLprob_plusplus = log(LLprob_plusplus);

                LL_plus_value = sum(LLprob_plus);
                LL_plusplus_value = sum(LLprob_plusplus);

                FisherIN(yf,yf) = (LL_plusplus_value-2*LL_plus_value+LL_value)/(precision*precision);

            end

            fprintf('Komentarz: Najpierw wyliczamy wartosci dla przekatnej\n')
            -FisherIN

            % pozostale wartosci

            ir=1;
            for ir=1:size(parameters_final,2)-1
                row_aux = ir;
                wg=1;
                for wg=1:size(parameters_final,2)
                    if wg>ir
                        parameters2_plusplus = parameters_final;
                        parameters2_plusminus = parameters_final;
                        parameters2_minusplus = parameters_final;
                        parameters2_minusminus = parameters_final;

                        parameters2_plusplus(:,row_aux) = parameters2_plusplus(:,row_aux)+precision;
                        parameters2_plusplus(:,wg) = parameters2_plusplus(:,wg)+precision;

                        parameters2_plusminus(:,row_aux) = parameters2_plusminus(:,row_aux)+precision;
                        parameters2_plusminus(:,wg) = parameters2_plusminus(:,wg)-precision;

                        parameters2_minusplus(:,row_aux) = parameters2_minusplus(:,row_aux)-precision;
                        parameters2_minusplus(:,wg) = parameters2_minusplus(:,wg)+precision;

                        parameters2_minusminus(:,row_aux) = parameters2_minusminus(:,row_aux)-precision;
                        parameters2_minusminus(:,wg) = parameters2_minusminus(:,wg)-precision;

                        V2_plusplus = attributes_final.*parameters2_plusplus;
                        V2_plusminus = attributes_final.*parameters2_plusminus;
                        V2_minusplus = attributes_final.*parameters2_minusplus;
                        V2_minusminus = attributes_final.*parameters2_minusminus;

                        tf = 1;
                        LL2_plusplus = [];
                        LL2_plusplus_max = [];
                        LL2_plusminus = [];
                        LL2_plusminus_max = [];
                        LL2_minusplus = [];
                        LL2_minusplus_max = [];
                        LL2_minusminus = [];
                        LL2_minusminus_max = [];
                        for tf=1:size(parameters_final,1)
                            aux_LL = tf;
                            se=1;
                            for se=1:nalt_aux
                                LL2_plusplus(aux_LL,se)=sum(V2_plusplus(aux_LL,alt_number_decoded2==se));
                                LL2_plusplus_max(aux_LL,1)=max(LL2_plusplus(aux_LL,:));

                                LL2_plusminus(aux_LL,se)=sum(V2_plusminus(aux_LL,alt_number_decoded2==se));
                                LL2_plusminus_max(aux_LL,1)=max(LL2_plusminus(aux_LL,:));

                                LL2_minusplus(aux_LL,se)=sum(V2_minusplus(aux_LL,alt_number_decoded2==se));
                                LL2_minusplus_max(aux_LL,1)=max(LL2_minusplus(aux_LL,:));

                                LL2_minusminus(aux_LL,se)=sum(V2_minusminus(aux_LL,alt_number_decoded2==se));
                                LL2_minusminus_max(aux_LL,1)=max(LL2_minusminus(aux_LL,:));
                            end
                        end

                        % finding max LL
                        gi=1;
                        choice2_plusplus_aux = [];
                        Choice2_plusplus = zeros(size(parameters_final,1),nalt_aux);;
                        choice2_plusminus_aux = [];
                        Choice2_plusminus = zeros(size(parameters_final,1),nalt_aux);;
                        choice2_minusplus_aux = [];
                        Choice2_minusplus = zeros(size(parameters_final,1),nalt_aux);;
                        choice2_minusminus_aux = [];
                        Choice2_minusminus = zeros(size(parameters_final,1),nalt_aux);;
                        for gi=1:size(parameters_final,1)
                            choice2_plusplus_aux(gi,1)=find(LL2_plusplus(gi,:)==LL2_plusplus_max(gi,1),1,'first');
                            Choice2_plusplus(gi,choice2_plusplus_aux(gi,1))=1;

                            choice2_plusminus_aux(gi,1)=find(LL2_plusminus(gi,:)==LL2_plusminus_max(gi,1),1,'first');
                            Choice2_plusminus(gi,choice2_plusminus_aux(gi,1))=1;

                            choice2_minusplus_aux(gi,1)=find(LL2_minusplus(gi,:)==LL2_minusplus_max(gi,1),1,'first');
                            Choice2_minusplus(gi,choice2_minusplus_aux(gi,1))=1;

                            choice2_minusminus_aux(gi,1)=find(LL2_minusminus(gi,:)==LL2_minusminus_max(gi,1),1,'first');
                            Choice2_minusminus(gi,choice2_minusminus_aux(gi,1))=1;
                        end

                        LL2_plusplus=exp(LL2_plusplus);
                        LL2_plusplus_up = LL2_plusplus.*Choice2_plusplus;

                        LL2_plusminus=exp(LL2_plusminus);
                        LL2_plusminus_up = LL2_plusminus.*Choice2_plusminus;

                        LL2_minusplus=exp(LL2_minusplus);
                        LL2_minusplus_up = LL2_minusplus.*Choice2_minusplus;

                        LL2_minusminus=exp(LL2_minusminus);
                        LL2_minusminus_up = LL2_minusminus.*Choice2_minusminus;

                        ct=1;
                        LL2prob_plusplus = [];
                        LL2prob_plusminus = [];
                        LL2prob_minusplus = [];
                        LL2prob_minusminus = [];
                        for ct= 1:size(parameters_final,1)
                            LL2prob_plusplus(ct,1) = sum(LL2_plusplus_up(ct,:))/sum(LL2_plusplus(ct,:));
                            LL2prob_plusminus(ct,1) = sum(LL2_plusminus_up(ct,:))/sum(LL2_plusminus(ct,:));
                            LL2prob_minusplus(ct,1) = sum(LL2_minusplus_up(ct,:))/sum(LL2_minusplus(ct,:));
                            LL2prob_minusminus(ct,1) = sum(LL2_minusminus_up(ct,:))/sum(LL2_minusminus(ct,:));
                        end

                        LL2prob_plusplus = log(LL2prob_plusplus);
                        LL2prob_plusminus = log(LL2prob_plusminus);
                        LL2prob_minusplus = log(LL2prob_minusplus);
                        LL2prob_minusminus = log(LL2prob_minusminus);

                        LL2_plusplus_value = sum(LL2prob_plusplus);
                        LL2_plusminus_value = sum(LL2prob_plusminus);
                        LL2_minusplus_value = sum(LL2prob_minusplus);
                        LL2_minusminus_value = sum(LL2prob_minusminus);

                        FisherIN(row_aux,wg) = (LL2_plusplus_value-LL2_plusminus_value-LL2_minusplus_value+LL2_minusminus_value)/(4*precision*precision);
                        FisherIN(wg,row_aux) = FisherIN(row_aux,wg);
                    end
                end
            end

    %         fractional_factorial_random
    %         attributes_final
    %         parameters_final

            fprintf('Komentarz: Oraz pozostale wartosci: \n')
            FisherIN=-FisherIN

            fprintf('Komentarz: AVC matrix to odwrocenie macierzy Fischera: inv() lub ^(-1)\n')
            AVCmatrix=inv(FisherIN) %^(-1)
            fprintf('Komentarz: Problemem jest to, ze czasem na przekatnej macierzy AVC wystepuja minusowe wartosci \n')
            diag_avc = diag(AVCmatrix)'
            fprintf('Komentarz: Nie powinno tak byc, gdyz standard_errors = sqrt(diag_avc)\n')
            standard_errors = sqrt(diag_avc)
            parameters_final;
            parameters_mean = sum(parameters_final,1) ./ sum(parameters_final~=0,1);

            fprintf('Komentarz: Nastepnie na podstawie macierzy AVC wyliczane sa odpowiednie miary\n')
            fprintf('Komentarz: W tym D-error, A-error, S-estimate i Utility Balance')
            % the lower the better
            Derror = det(AVCmatrix)^(1/size(parameters_final,2))
            Aerror = trace(AVCmatrix)/size(parameters_final,2)
            Sestimate_all = ((1.96 .* standard_errors)./parameters_mean).^2
            Sestimate = max(Sestimate_all)

            % utility balance
            for p= 1:size(parameters_final,1)
                aux_ub = p;
                for o=1:nalt_aux
                    LL_ub(p,o) = LL(p,o)/sum(LL(p,:)); 
                end
            end

            for b=1:size(parameters_final,1)
                aux2_ub = b;
                for k=1:nalt_aux
                    LL_ub2(b,k)=LL_ub(b,k)/(1/nalt_aux);
                end
            end

            LL_ub3 = prod(LL_ub2,2);
            UB = mean(LL_ub3)

            fprintf('Komentarz: W zaleznie od zadeklarowanego kryterium funkcja szuka najmniejszej wartosci (ew. najwiekszej dla B)\n')
            fprintf('Komentarz: Mozna zadeklarowac ilosc powtorzen, funkcja w kazdym powtorzeniu robi te same wyliczenia\n')
            fprintf('Komentarz: Jesli znajduje lepszy design to podmienia go w pamieci\n')
            if min_cond==["A"]
                    if Aerror<Aerror_aux
                        Aerror_aux=Aerror;
                        Derror_aux=Derror;
                        UB_aux = UB;
                        Sestimate_all_aux = Sestimate_all;
                        Sestimate_aux = Sestimate;
                        fractional_factorial_random_aux = fractional_factorial_random;
                        fprintf('Iteration %d:   A-error: %f,   D-error: %f,   Utility Balance: %f,   S-estimate: %f\n', h, Aerror, Derror, round(UB,4), Sestimate)
                    end
            elseif min_cond==["D"]
                    if Derror<Derror_aux
                        Aerror_aux=Aerror;
                        Derror_aux=Derror;
                        UB_aux = UB;
                        Sestimate_all_aux = Sestimate_all;
                        Sestimate_aux = Sestimate;
                        fractional_factorial_random_aux = fractional_factorial_random;
                        fprintf('Iteration %d:   A-error: %f,   D-error: %f,   Utility Balance: %f,   S-estimate: %f\n', h, Aerror, Derror, round(UB,4), Sestimate)
                    end
            elseif min_cond==["B"]
                    if UB>UB_aux
                        Aerror_aux=Aerror;
                        Derror_aux=Derror;
                        UB_aux = UB;
                        Sestimate_all_aux = Sestimate_all;
                        Sestimate_aux = Sestimate;
                        fractional_factorial_random_aux = fractional_factorial_random;
                        fprintf('Iteration %d:   A-error: %f,   D-error: %f,   Utility Balance: %f,   S-estimate: %f\n', h, Aerror, Derror, round(UB,4), Sestimate)
                    end
            elseif min_cond==["S"]
                    if Sestimate<Sestimate_aux
                        Aerror_aux=Aerror;
                        Derror_aux=Derror;
                        UB_aux = UB;
                        Sestimate_all_aux = Sestimate_all;
                        Sestimate_aux = Sestimate;
                        fractional_factorial_random_aux = fractional_factorial_random;
                        fprintf('Iteration %d:   A-error: %f,   D-error: %f,   Utility Balance: %f,   S-estimate: %f\n', h, Aerror, Derror, round(UB,4), Sestimate)
                    end                
            else
                cprintf(rgb('Red'),'Error: Incorrect argument in Design.EF_errors! (set it to ["A"] for A-error, ["D"] for D-error), ["B"] for Utility Balance or ["S"] for S-estimate \n')
            end
            % singular matrix
       end 
   end
   
   

   for iopt=1:size(vector_name,2)
       fractional_factorial_random_aux(fractional_factorial_random_aux==strcat(num2str(iopt),'*'))=vector_name(iopt);
   end
   fractional_factorial_random_aux;
   fprintf('Komentarz: Na sam koniec znamy najlepszy znaleziony dotychczas design')
   fractional_factorial_random_aux
end