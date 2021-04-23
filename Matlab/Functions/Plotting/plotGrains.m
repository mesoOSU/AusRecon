function plotGrains(Grains,label)

figure; plot(Grains.grains)
if label == 1
    text(Grains.grains,num2str(Grains.grains.id))
end
end