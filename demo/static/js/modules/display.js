async  function updateInput(stepId) {
    steps[stepId].params = document.querySelector(`#input-step-${stepId} input`).value;
    updateOutputArea();
}

let isUpdating = false;
async function updateOutputArea() {
    if (isUpdating) { console.log('Update is currently in progress. Please wait.'); return;}
    isUpdating = true;
    const outputArea = document.getElementById('outputArea');
    outputArea.classList.add('custom-border'); // 应用自定义边框样式
    outputArea.innerHTML = '<h2>Loading...</h2>';
    for (const step of steps) {
        try {
            step.params = document.querySelector(`#input-step-${step.id} input`).value;
            const response = await fetch('/common/get_display_value', {
                method: 'POST',  headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ mode: step.mode, stepId: step.id,  params: step.params, step_type: step.type})
            });
            const data = await response.json();
            const outputDiv = document.createElement('div');
            outputDiv.id = `output-step-${step.id}`;
            outputDiv.classList.add('alert', 'custom-alert-info'); // 应用自定义警告信息样式
            outputDiv.innerHTML = `<strong>${step.type}:</strong> ${data.displayString}`;
            outputArea.appendChild(outputDiv);
        } catch (error) { console.error('Error:', error); outputArea.innerHTML = '<h2>Error Loading Data</h2>'; }
    }
    isUpdating = false;
}


