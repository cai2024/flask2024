

async function initialPreset(presetType) {
    try {
        const response = await fetch('/get_preset', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ presetType: presetType })
        });
        const result = await response.json();
	const preset_options = result.options;
	const preset_select = result.select;
	const preset_input = result.params;

	clickButtons(preset_options).then(() => {
        presetSteps(preset_select,preset_input);
        })
	disableAllButtons()
	


    } catch (error) {
        console.error('Error initialPreset:', error);
    }
}


function clickButtons(soft) {
    let promises = soft.map((element, index) => {
        return new Promise((resolve) => {
            setTimeout(() => {
                addStep(element);
                resolve(); 
            }, index * 10);
        });
    });
    return Promise.all(promises);
}



// 初始化
async function presetSteps(preset_select,preset_input) {
    confirmBtn.disabled = true;  // 禁用按钮
    for (let i = 0; i < steps.length; i++) {
        steps[i].id = i;
    }
    for (const step of steps) {
        await presetOptions(step.type, step.id,preset_select,preset_input);
    }
    $('#selected-steps').sortable('disable');
    const submitBtn = document.getElementById('submitBtn');
    submitBtn.disabled = false; 
}




async function presetOptions(stepType, stepId,preset_select,preset_input) {
    try {
        const response = await fetch('/get_options', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ step_type: stepType, my_options: myOptions })
        });
        const data = await response.json();
        const options = data.options;  // 核心任务，传递给更新UI
        await presetOptionsUI(stepId, stepType, options, preset_select,preset_input);   //核心修改
        // 检查 stepType 的前三个字符是否为 "map"
        if (stepType.slice(0, 3) === 'map') {
            const mapResponse = await fetch('/index_select', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ step_type: stepType })
            });
            const mapData = await mapResponse.json();
            const index_options = mapData.index_options;
            await initIndex(stepId, stepType, index_options);
        }
    } catch (error) {
        console.error('Error:', error);
    }
}
async function presetOptionsUI(stepId, stepType, options, preset_select,preset_input) {
    const step = steps.find(step => step.id === stepId);
    for (let i = 0; i < options.length; i++) {
	step.params = preset_input[stepId] 
    }
    // 核心修改，通过列表索引假定对应关系，这也意味着先点击按钮会报错，可以解决，先传入当前index就行了
    const optionsHtml = options.map(option => `<option value="${option}">${option}</option>`).join('');
    const inputDiv = document.createElement('div');
    inputDiv.id = `input-step-${stepId}`;
    inputDiv.classList.add('form-group');
    inputDiv.innerHTML = `
        <label><strong>step${stepId}: ${stepType}:</strong></label>
        <input id="input-${stepId}"  type="text" class="form-control" placeholder="请输入内容" onchange="updateInput(${stepId})" value="${steps.find(step => step.id === stepId).params || ''}">
        <select id="select-${stepId}" class="form-control mt-2" onchange="updateSelectOption(${stepId})">${optionsHtml}</select>
    `;
    document.getElementById('inputArea').appendChild(inputDiv);
    // 绑定输入框
    const inputElement = document.querySelector(`#input-step-${stepId} input`);
    step.params = inputElement.value;
    inputElement.addEventListener('input', () => {
        updateOutputArea();
    });
    // 绑定选择框
    const selectElement = document.querySelector(`#input-step-${stepId} select`);
    if (selectElement && options.length > 0) {
        for (let i = 0; i < options.length; i++) {        
        if (options[i] === preset_select[stepId]) {
            selectElement.selectedIndex = i;
	    break;
	} }
        step.select_option = selectElement.value;
	// 核心修改，遍历当前选项框获取，列表index的对应位置，问题同上
        await updateMyOptions(stepType, stepId);
    }else {
    console.log("选项匹配失败");
}
}


