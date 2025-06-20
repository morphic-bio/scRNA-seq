#include "queue.h"
#include <stdio.h>

// Function to initialize the queue
void init_queue(Queue *queue) {
    queue->capacity = INITIAL_CAPACITY;
    queue->data = (uint64_t *)malloc(queue->capacity * sizeof(uint64_t));
    queue->front = 0;
    queue->back = 0;
    queue->size = 0;
}

void clear_queue(Queue *queue) {
    queue->front = 0;
    queue->back = 0;
    queue->size = 0;
}
// Function to check if the queue is empty
int is_empty(Queue *queue) {
    return queue->size == 0;
}

// Function to get the number of elements in the queue
size_t size(Queue *queue) {
    return queue->size;
}

// Function to enqueue an element to the back of the queue
void enqueue(Queue *queue, uint64_t value) {
    // Check if the queue needs to be resized
    if (queue->size == queue->capacity) {
        queue->capacity *= 2;
        queue->data = (uint64_t *)realloc(queue->data, queue->capacity * sizeof(int));

        // Adjust the front and back indices if the queue is wrapped around
        if (queue->front > queue->back) {
            for (int i = 0; i < queue->front; i++) {
                queue->data[queue->capacity / 2 + i] = queue->data[i];
            }
            queue->back = queue->capacity / 2 + queue->front;
        }
    }

    // Enqueue the element
    queue->data[queue->back] = value;
    queue->back = (queue->back + 1) % queue->capacity;
    queue->size++;
}

// Function to dequeue an element from the front of the queue
uint64_t dequeue(Queue *queue) {
    if (is_empty(queue)) {
        printf("Queue underflow\n");
        exit(EXIT_FAILURE);
    }

    uint64_t value = queue->data[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size--;
    return value;
}

// Function to get the element at the front of the queue without dequeuing it
uint64_t peek(Queue *queue) {
    if (is_empty(queue)) {
        printf("Queue is empty\n");
        exit(EXIT_FAILURE);
    }
    return queue->data[queue->front];
}

// Function to free the memory allocated for the queue
void free_queue(Queue *queue) {
    free(queue->data);
    queue->data = NULL;
    queue->front = 0;
    queue->back = 0;
    queue->size = 0;
    queue->capacity = 0;
}
