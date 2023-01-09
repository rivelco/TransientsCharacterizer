clear, clc

disp("Working...")

load("crf\Stoixeion_01_04_CRFS_K3 UDF 10\results.mat")
load("crf\Stoixeion_01_04_CRFS_K3 UDF 10\best_model.mat")
load("dbs\Stoixeion_01_04_CRFS_K3.mat")

%import relevant information
num_neur = size(data,2);
num_ens = size(UDF, 2);
ens = results.core_crf;
ENS_STATE = cell(1,num_ens);

node_str=results.epsum(1:num_neur);
degrees = sum(best_model.structure(1:num_neur,1:num_neur));
auc=results.auc(1:num_neur,:);
edge_potentials = best_model.theta.edge_potentials(1:num_neur,1:num_neur);

for i = 1:num_ens
    ENS_STATE{i}=auc(:,i);
    ENS_STATE{i}=round(ENS_STATE{i});
    size_ens{i}=sum(ENS_STATE{i});
    ENS_STATE{i}=ENS_STATE{i}.*transpose((1:num_neur));
    ENS_STATE{i}(ENS_STATE{i}==0)=[];
end

nsmi = []; %min node str
nsma = []; %max node str
nS = []; %min max norm node str
auc2 = []; %auc by ens

PAPS_INDEXED = cell(2,num_ens); %INDEX OF PAPS BY ENS
nodes_strength = cell(1, num_ens);
AUCs = cell(1, num_ens);
for i = 1:num_ens
    
    nsmi= min(node_str([ENS_STATE{i}]));
    nsma = max(node_str([ENS_STATE{i}]));
    nS = transpose((node_str([ENS_STATE{i}])-nsmi)./(nsma-nsmi));

    auc2 = auc([ENS_STATE{i}],i);
    auc2=transpose(auc2);
    
    PAPS_INDEXED{1,i}=ENS_STATE{i};
    for ii = 1:length(ENS_STATE{i})
        PAPS = @(Ni) (nS(Ni)+auc2(Ni))/2;
        PAPS_INDEXED{2,i}=[PAPS_INDEXED{2,i} PAPS(ii)];
        PAPS_INDEXED{2,i}(isnan(PAPS_INDEXED{2,i}))=0;
    end

    cellInEnsamble = transpose(ENS_STATE{i});
    %nS(isnan(nS)) = 0;
    nodes_strength{i} = [cellInEnsamble; nS];

    %auc2(isnan(auc2)) = 0;
    AUCs{i} = [cellInEnsamble; auc2];
end

save("CRF_Stoixeion_01_04_UDF_phi10.mat", "nodes_strength", "AUCs")
disp("Done.")
