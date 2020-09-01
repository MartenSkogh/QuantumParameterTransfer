
class Callbacker():

    def __init__(self, filename, hf_energy=0):
        self.output_file = open(filename, 'w')
        self.hf_energy = hf_emergy


    def callback(self, eval_count, param_set, mean, estimator_error):
        # You can overwrite this if you want to, and have another function
        string_rep = ','.join([str(eval_count), str(self.hf_energy + mean)] + [str(i) for i in param_set])
        #self.output_file.write(f'{string_rep}\n')
        print(string_rep)

    def close(self):
        self.output_file.close()


