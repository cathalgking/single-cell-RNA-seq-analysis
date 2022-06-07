def label2index(label):
    """ 
    Return the appropriately formatted cell_index based on input cell_label
    
    """
    if '-' in label:
        cell_label = [int(n) for n in label.split('-')]
        return str((cell_label[0] - 1) * 96 * 96 + (cell_label[1] - 1) * 96 + cell_label[2])
    else:
        return str(label)


print(label2index('1-1-1'))
print(label2index('5-5-5'))
print(label2index('43-12-77'))
print(label2index('96-96-96'))
print(label2index('5-43-9'))