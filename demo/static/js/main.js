// 公用的更新全局my_options
async function updateMyOptions(stepType, stepId) {
    const step = steps.find(step =>  step.id === stepId);
    try {
        const response = await fetch('/update_my_options', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ step_id: stepId, step_type: stepType, select_option: step.select_option, my_options: myOptions })
        });
        const data = await response.json();
        const updatedMYOptions = data.my_options;  // 实现初始化时更新my_options
        myOptions = updatedMYOptions;
    } catch (error) {
        console.error('Error:', error);
    }
}


let isUpdating = false;
async function updateOutputArea() {
    if (isUpdating) {
        console.log('Update is currently in progress. Please wait.');
        return;
    }
    isUpdating = true;
    const outputArea = document.getElementById('outputArea');
    outputArea.classList.add('custom-border'); // 应用自定义边框样式
    outputArea.innerHTML = '<h2>Loading...</h2>';

    for (const step of steps) {
        try {
            const inputElement = document.querySelector(`#input-step-${step.id} input`);
            const inputValue = inputElement.value;
            step.params = inputValue;
            const response = await fetch('/get_display_value', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    type: step.type,
                    select_option: step.select_option,
                    params: step.params || ""
                })
            });
            const data = await response.json();
            const outputDiv = document.createElement('div');
            outputDiv.id = `output-step-${step.id}`;
            outputDiv.classList.add('alert', 'custom-alert-info'); // 应用自定义警告信息样式
            outputDiv.innerHTML = `<strong>${step.type}:</strong> ${data.displayString}`;
            outputArea.appendChild(outputDiv);
        } catch (error) {
            console.error('Error:', error);
            outputArea.innerHTML = '<h2>Error Loading Data</h2>';
        }
    }
    isUpdating = false;
}
