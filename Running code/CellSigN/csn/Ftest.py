import numpy as np

# Define two lists
list1 = np.array([2, 3, 5, 6, 8, 9, 11, 12, 14, 15])
list2 = np.random.randint(1, 100, size=90)

# Calculate individual variances
mean_list1 = np.mean(list1)
mean_list2 = np.mean(list2)
variance_list1 = np.var(list1, ddof=1)
variance_list2 = np.var(list2, ddof=1)
n1 = len(list1)
n2 = len(list2)

# Calculate combined variance
combined_variance = (((n1 - 1) * variance_list1 + (n2 - 1) * variance_list2 + (n1 * n2) / (n1 + n2)
                      * (mean_list1 - mean_list2)**2) / (n1 + n2 - 1))

# Output the results
print("Variance of list1:", variance_list1)
print("Variance of list2:", variance_list2)
print("Combined variance:", combined_variance)