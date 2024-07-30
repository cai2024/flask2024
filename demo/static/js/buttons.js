let steps = [];

function addStep(type) {
    const stepId = steps.length;
    const stepItem = `<li class="list-group-item step-item" id="step-${stepId}">${type}</li>`;
    $('#selected-steps').append(stepItem);
    $('#selected-steps').sortable({
        update: function(event, ui) {
            updateStepOrder();
        }
    });
    $('#selected-steps').disableSelection();
    steps.push({ id: stepId, type: type, text: '', select_option: '' });
}

function updateStepOrder() {
    const newOrder = [];
    $('#selected-steps').children().each(function(index, element) {
        const stepId = parseInt(element.id.replace('step-', ''));
        newOrder.push(steps.find(step => step.id === stepId));
    });
    steps = newOrder;
}

