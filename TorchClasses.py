from torch import nn

class BiologicalSequenceDataset:
    def __init__(self, sequenceData, rxnData):
        self.seqData = sequenceData
        self.rxnData = rxnData
    
    def __len__(self):
        return len(self.seqData)
    
    def __getitem__(self, i):
        seq = self.seqData[i]
        rxn = self.rxnData[i]
        return seq, rxn

class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(2048, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 7182),
            nn.ReLU(),
            nn.Linear(7182, 7182)
        )
    
    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return(logits)

