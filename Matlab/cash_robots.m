close all
clear all

days = 380;

energy_per_hour = [7 72 375 1920 9600 25000 60000 120000];
cost_per_bot = [100 1000 5000 25000 100000 250000 500000 900000 ]; % in diamonds
owned_bots = zeros(days,8);
owned_bots(1,:) = [1 0 0 0 0 0 0 0];
energy_diamond_ratio = 100;
bonus_rate=50;

%plot((energy_per_hour./cost_per_bot));

% First Scenario - Always Buy first Robot when Possible
energy = 0;
diamonds = zeros(days,1);
wallet=zeros(days,1);
for j=1:5%length(cost_per_bot)
bot_to_buy = j;
for i=1:days
    energy = energy + sum(owned_bots(i,:).*energy_per_hour)*24;
    if energy >= 10000
        energy_sold = energy-mod(energy,300);
        energy = mod(energy,300);
        converted_diamonds = energy_sold/energy_diamond_ratio; 
        earned_diamonds = floor(converted_diamonds*.7);
        wallet_diamonds = converted_diamonds-earned_diamonds;
        earned_diamonds = earned_diamonds +floor(rand()*bonus_rate);
        %fprintf('%d\n',earned_diamonds);
    else
        wallet_diamonds = 0;
        earned_diamonds = floor(rand()*bonus_rate);
    end
    wallet(i+1) = wallet(i)+wallet_diamonds;
    diamonds(i)=diamonds(i)+earned_diamonds;
    
    if diamonds(i) >= cost_per_bot(bot_to_buy)
        owned_bots(i+1,bot_to_buy) = owned_bots(i,bot_to_buy) + (diamonds(i)-mod(diamonds(i),cost_per_bot(bot_to_buy)))/cost_per_bot(bot_to_buy);
        diamonds(i) = mod(diamonds(i),cost_per_bot(bot_to_buy));
        diamonds(i+1) = diamonds(i);
    else
        owned_bots(i+1,bot_to_buy)=owned_bots(i,bot_to_buy);
        diamonds(i+1) = diamonds(i);
    end
    
end
d(:,j)=diamonds;
owned(:,j)=owned_bots(:,bot_to_buy);
w(:,j)=wallet(:);
end
figure;
subplot(2,2,1)
plot(log10(d))
title('Diamonds')
subplot(2,2,2);
plot(log10(owned))
title('Owned Bots Type 1')

subplot(2,2,3);
plot(log10(w))
title('Diamonds to Sell');

owned(95,:)
    
    