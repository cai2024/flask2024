document.getElementById('submitBtn').addEventListener('click', submitSteps);
document.getElementById('unlockBtn').addEventListener('click', function() {initialPreset("WGS");});




// 上传fq_type, raw_myOptions[0], steps到session




async function submitSteps() {
    try {
        soft_step.forEach((step, stepId) => {
            step.params = document.getElementById(`input-${stepId}`).value

            const selectElements = document.querySelectorAll(`#input-step-${stepId} select`);
            selectElements.forEach(select => {
                step.select_option[select.title] = select.value;
            });

            const indexElements = document.querySelectorAll(`#index-step-${stepId} select`);
            indexElements.forEach(select => {
                step.select_index[select.title] = select.value;
            });

        })
        const response = await fetch('/batch_code/store_soft_step', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({soft_step: soft_step})
        });
        if (response.ok) {window.location.href = '/batch_code/submit_form_code/' + task_id;}

    }catch (error) {console.error('Error submitting steps:', error);}
}





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

