import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, EsmForMaskedLM
import sys
import csv

fasta_file = sys.argv[1] #the path of the single sequence fasta file
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 
               'I', 'K', 'L', 'M', 'N', 'P', 'Q', 
               'R', 'S', 'T', 'V', 'W', 'Y']

with open(fasta_file, 'r') as file:
    lines = file.readlines()
    protein_sequence = ''.join(line.strip() for line in lines[1:]) 

tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
model = EsmForMaskedLM.from_pretrained("facebook/esm2_t33_650M_UR50D")

#inputs = tokenizer("LIRD<mask>ISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE", return_tensors="pt")

input_ids = tokenizer.encode(protein_sequence, return_tensors="pt")
sequence_length = len(protein_sequence) #input_ids.shape[1] - 2

end_pos = sequence_length -1
start_pos = 0

all_probabilities = []
for position in range(start_pos, end_pos + 1):
    masked_input_ids = input_ids.clone()
    masked_input_ids[0, position] = tokenizer.mask_token_id
    with torch.no_grad():
        # Get logits from the model
        logits = model(masked_input_ids).logits

    # Apply softmax to convert logits to probabilities
    probabilities = F.softmax(logits[0, position], dim=0)
    all_probabilities.append(probabilities)

csv_fname = "ESM_probabilities_" + fasta_file + ".csv"

with open(csv_fname, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        header = ["Position"] + amino_acids
        csv_writer.writerow(header)
        for pos, probs in enumerate(probabilities):
            csv_writer.writerow([pos, *probs])
