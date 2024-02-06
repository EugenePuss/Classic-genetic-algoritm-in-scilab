//функция приспособленности
function adapted = f_adaptet(popul)
    
    adapted = zeros(size(popul));
    for i =1:size(popul,1)
        adapted(i) = sin(popul(i))*popul(i);
    end
endfunction

//функция селекции хромосом
function selection = f_select(pop,adapt_mass,num_sel)
    selection = zeros(num_sel,1);
    prob = adapt_mass/sum(adapt_mass);
    // создаем вектор суммированных вероятностей (для рулетки)
    cumulative_prob = cumsum(prob);
    for i = 1:num_sel
        r = rand();
        j = find(cumulative_prob >= r, 1);
        selection(i) = pop(j) ;
    end
endfunction

//функция скрещивания 
function offspring = crossover(parents, crossover_rate)
    offspring = zeros(size(parents, 1),1);
    for i = 1:2:size(parents,1)
        //выбираем родителей
         rand_1 = floor(rand()*size(parents,1))+1;
            par1 = parents(rand_1);
            parents(rand_1)=[];
            rand_2 = floor(rand()*size(parents,1))+1;
            par2 = parents(rand_2);
            parents(rand_2)=[];
        if rand() < crossover_rate
            //получаем потомков
            point_cross=grand(1,1,"uin",1,Chrom_s);
            offspring(i) = bitset(par1,[1:abs(point_cross- Chrom_s)],[bitget(par2,1:abs(point_cross- Chrom_s))]);
            offspring(i+1) = bitset(par2,[1:abs(point_cross- Chrom_s)],[bitget(par1,1:abs(point_cross- Chrom_s))]);
        else 
            offspring(i)= par1;
            offspring(i+1) = par2;
        end
    end
end

//функции мутации хромосом
function mutation = f_mut(popul,ver_mutation)
    mutation = zeros(size(popul,1),1);
    for i=1:size(popul,1) 
        if rand() <= ver_mutation
            point_mut=grand(1,1,"uin",1,Chrom_s);
            gen_mut = 0;
            if bitget(popul(i),abs(point_mut-Chrom_s)+1:abs((point_mut-Chrom_s)+1))==1
                get_mut = 0;
            else get_mut = 1
            end
            mutation(i)=bitset(popul(i),[point_mut],[gen_mut]);
        else 
            mutation(i)=popul(i);
        end
    end
endfunction

//Функция отбора наилучшей хромосомы
function best_ch = f_best(pop)
    best_fit = -50000;
    best_ch = zeros(size(pop));
    fit = f_adaptet(pop)
    disp(fit)
    for i =1:size(pop,1)
        disp(pop(i))
     if fit(i) > best_fit
         best_fit = fit(i);
         best_ch = pop(i);
     end 
    end
endfunction
// генетический алгоритм
PopSize=160;//// размер популяции
Crossing=0.9 //вероятность скрещивания
Mutation=0.5;//вероятность мутации
Pokol=5;//количество популяций 
Chrom_s = 0; //размер хромосомы
best=0;// самый приспособленный
// Инициализация начальной популяции
    population = grand(PopSize, 1, "uin", 1, 265);
    Chrom_s = length(dec2bin(max(population)));
    x=linspace(1,PopSize,100)// график функции 
    plot(x, f_adaptet(x'))
    xgrid(5), gca.children.children.thickness=3;
    for i = 1:Pokol
        // Оценка приспособленности каждого индивидуума
        fitness = f_adaptet(population);
        disp(fitness)
        // Выбор родителей для скрещивания /селекция
            parent = f_select(population,fitness,floor(size(fitness,2)/2));
            disp(parent)
             disp(dec2bin(parent))
        // Скрещивание родителей
        offsp = crossover(parent, Crossing);
        disp(offsp)
         //
         disp(dec2bin(offsp))
        // Мутация потомства
        mutated_offspring = f_mut(offsp, Mutation);
        disp( mutated_offspring)
        population = mutated_offspring;
    end
    // Нахождение лучшего индивидуума в популяции
    best = f_best(population)
    disp(best)
    disp(dec2bin( best))
    
    

