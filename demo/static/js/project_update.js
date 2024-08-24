async  function updateInput(stepId) {
    const inputElement = document.querySelector(`#input-step-${stepId} input`);
    const inputValue = inputElement.value;
    steps[stepId].params = inputValue;
    updateOutputArea();
}


// 更新选项
async  function updateSelectOption(stepId) {
    const step = steps.find(step => step.id === stepId);
    const key_dir={}
    for (let find_key in step.select_option) {
        const selectElement = document.querySelector(`#input-step-${stepId} select#select-${stepId}-${find_key}`)
        if (step.select_option[find_key] !== selectElement.value){
            const change_key= find_key
            step.select_option[find_key] = selectElement.value
            try {
                for (let key in step.select_option) {
                    const selectElement = document.querySelector(`#input-step-${stepId} select#select-${stepId}-${key}`)
                    const optionsList = Array.from(selectElement.options).map(option => option.text);
                    key_dir[key]=[step.select_option[key],optionsList]}
                const response = await fetch('/update_pair_select', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({key_dir: key_dir,change_key: change_key,step_type:steps.type,mode:step.mode})
                });
                const update_select_dir = (await response.json()).update_select_dir;
                for (let key in update_select_dir) {
                    const selectElement = document.querySelector(`#input-step-${stepId} select#select-${stepId}-${key}`)
                    const optionsList = Array.from(selectElement.options).map(option => option.text);
                    const index = optionsList.indexOf(update_select_dir[key][0]);
                    if (index !== -1) {selectElement.selectedIndex = index;}
                    step.select_option[key] = selectElement.value
                }
            } catch (error) { console.error('Error:', error); }
            break
        }
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

