# commonly used mpi functions
import numpy

class MPIRoutines:

    @staticmethod
    def gather_single_value(value, number_of_values, comm, my_rank):

        temp_full_value = comm.gather(value, root=0)

        if my_rank == 0:
            full_value = numpy.zeros(number_of_values, dtype=numpy.float64)

            # putting all the obtained minimum values into a single array
            count = 0
            for i in temp_full_value:
                for j in i:
                    full_value[count] = j
                    count += 1
        else:
            full_value = None

        return full_value

    @staticmethod
    def gather_list_of_values(values, comm, my_rank):

        full_values = comm.gather(values, root=0)

        if my_rank == 0:
            list_of_values = []
            for i in range(len(full_values)):
                list_of_values += full_values[i]
        else:
            list_of_values = []

        return list_of_values

    @staticmethod
    def gather_numpy_array_of_values(values, comm, my_rank):

        full_values = comm.gather(values, root=0)

        if my_rank == 0:
            list_of_values = []
            for i in range(len(full_values)):
                list_of_values += full_values[i]
            array_of_values = numpy.zeros((len(list_of_values), list_of_values[0].shape[0]), dtype=numpy.float64)
            for i in range(len(list_of_values)):
                array_of_values[i, :] = list_of_values[i]
        else:
            array_of_values = None

        return array_of_values

    @staticmethod
    def distribute_points(samples, my_rank, num_cores, comm):

        sample_portion = None

        if my_rank == 0:

            # number of tasks per core
            tasks = len(samples) // num_cores  # // calculates the floor

            # remainder
            r = len(samples) - num_cores * tasks

            # array that holds how many tasks each core has
            tasks_core = numpy.zeros(num_cores, dtype=numpy.int64)
            tasks_core.fill(tasks)

            # distributing in the remainder
            ii = 0
            while r > 0:
                tasks_core[ii] += 1
                r -= 1
                ii += 1

            sample_portion = samples[0:tasks_core[0], :]

            if num_cores > 1:
                for i in range(1, num_cores):
                    start = sum(tasks_core[0:i])
                    end = start + tasks_core[i]
                    comm.send(samples[start:end, :], dest=i, tag=i * 11)

        else:
            if num_cores > 1:
                sample_portion = comm.recv(source=0, tag=my_rank * 11)

        return sample_portion

    @staticmethod
    def distribute_list_of_points(samples, my_rank, num_cores, comm):

        sample_portion = None

        if my_rank == 0:

            # number of tasks per core
            tasks = len(samples) // num_cores  # // calculates the floor

            # remainder
            r = len(samples) - num_cores * tasks

            # array that holds how many tasks each core has
            tasks_core = numpy.zeros(num_cores, dtype=numpy.int64)
            tasks_core.fill(tasks)

            # distributing in the remainder
            ii = 0
            while r > 0:
                tasks_core[ii] += 1
                r -= 1
                ii += 1

            sample_portion = samples[0:tasks_core[0]]

            if num_cores > 1:
                for i in range(1, num_cores):
                    start = sum(tasks_core[0:i])
                    end = start + tasks_core[i]
                    comm.send(samples[start:end], dest=i, tag=i * 11)

        else:
            if num_cores > 1:
                sample_portion = comm.recv(source=0, tag=my_rank * 11)

        return sample_portion