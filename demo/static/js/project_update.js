async  function updateInput(stepId) {
    const inputElement = document.querySelector(`#input-step-${stepId} input`);
    const inputValue = inputElement.value;
    steps[stepId].params = inputValue;
    updateOutputArea();
}


// 更新选项
async  function updateSelectOption(stepId) {
    const selectElement = document.querySelector(`#input-step-${stepId} select`);
    const selectOption = selectElement.value;
    steps[stepId].select_option = selectOption;
    myOptions=[raw_myOptions[0]];  // 后面要修改为初始化状态
    for (let i = 0; i < stepId + 1; i++) {
        await updateMyOptions(steps[i].type, steps[i].id); //调用公用更新函数
    }
    for (let i = stepId + 1; i < steps.length; i++) {
        await updateOptions_select(steps[i].type, steps[i].id);
    }
}
async function updateOptions_select(stepType, stepId) {
    try {
        const step = steps.find(step => step.id === stepId);
        const response = await fetch('/get_options', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ step_type: stepType, my_options: myOptions })
        });
        const data = await response.json();
        const options = data.options;

        const optionsHtml = options.map(option => `<option value="${option}">${option}</option>`).join('');
        const selectElement = document.querySelector(`#input-step-${stepId} select`);
        selectElement.innerHTML = optionsHtml;
        if (selectElement && options.length > 0) {
            const longestIndex = options.reduce((maxIndex, option, index, array) => option.length > array[maxIndex].length ? index : maxIndex, 0);
        selectElement.selectedIndex = longestIndex;
            step.select_option = selectElement.value; // 更新当前行选项，这里应该和初始化的保持一致
            await updateMyOptions(stepType, stepId); //提交给后台更新my_options
        }
    } catch (error) {
        console.error('Error:', error);
    }
}












// 展示部分的更新

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

